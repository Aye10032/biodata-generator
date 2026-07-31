import shutil
import tempfile
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from bioflow_sim.core.config import BatchCase, BatchConfig
from bioflow_sim.simulators.bulk_rna import simulate_bulk_rna
from bioflow_sim.simulators.chromatin import simulate_chromatin
from bioflow_sim.simulators.dna import simulate_dna
from bioflow_sim.simulators.hic import simulate_hic
from bioflow_sim.simulators.methylation import simulate_methylation
from bioflow_sim.simulators.scrna import simulate_scrna
from bioflow_sim.simulators.tumor_normal import simulate_tumor_normal

Simulator = Callable[..., dict[str, object]]

SIMULATORS: dict[str, Simulator] = {
    'bulk-rna': simulate_bulk_rna,
    'chromatin': simulate_chromatin,
    'dna': simulate_dna,
    'hic': simulate_hic,
    'methylation': simulate_methylation,
    'scrna': simulate_scrna,
    'tumor-normal': simulate_tumor_normal,
}


@dataclass(frozen=True)
class BatchResult:
    case_name: str
    simulator: str
    status: str
    message: str


def select_cases(
    config: BatchConfig,
    selected_names: Iterable[str],
) -> tuple[BatchCase, ...]:
    requested = tuple(selected_names)
    available = {case.name: case for case in config.cases}
    unknown = sorted(set(requested) - set(available))
    if unknown:
        raise ValueError(f'unknown case name(s): {", ".join(unknown)}')
    if requested:
        return tuple(available[name] for name in requested)
    return tuple(case for case in config.cases if case.enabled)


def run_batch(
    config: BatchConfig,
    *,
    selected_names: Iterable[str] = (),
    dry_run: bool = False,
    continue_on_error: bool = False,
) -> list[BatchResult]:
    cases = select_cases(config, selected_names)
    results: list[BatchResult] = []
    for case in cases:
        simulator = SIMULATORS.get(case.simulator)
        if simulator is None:
            error = f'unknown simulator {case.simulator!r}'
            if not continue_on_error:
                raise ValueError(f'{case.name}: {error}')
            logger.error('{}: {}', case.name, error)
            results.append(BatchResult(case.name, case.simulator, 'failed', error))
            continue

        output_dir = case.parameters['output_dir']
        if dry_run:
            message = f'would publish {config.publish} output to {output_dir}'
            logger.info(
                '[dry-run] {} ({}, {}) -> {}',
                case.name,
                case.simulator,
                config.publish,
                output_dir,
            )
            results.append(BatchResult(case.name, case.simulator, 'planned', message))
            continue

        logger.info('Starting case {} using {}', case.name, case.simulator)
        try:
            _run_case(simulator, case.parameters, config.publish)
        except (OSError, TypeError, ValueError) as error:
            if not continue_on_error:
                raise ValueError(f'{case.name}: {error}') from error
            logger.exception('Case {} failed: {}', case.name, error)
            results.append(BatchResult(case.name, case.simulator, 'failed', str(error)))
        else:
            results.append(BatchResult(case.name, case.simulator, 'completed', str(output_dir)))
    return results


def _run_case(
    simulator: Simulator,
    parameters: dict[str, object],
    publish: str,
) -> None:
    if publish == 'development':
        simulator(**parameters)
        return

    output_dir = parameters['output_dir']
    if not isinstance(output_dir, Path):
        raise TypeError('output_dir must be a path')
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise FileExistsError(f'output directory is not empty: {output_dir}')
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=f'.{output_dir.name}.bioflow-sim-',
        dir=output_dir.parent,
    ) as temporary:
        staging_dir = Path(temporary) / 'case'
        staging_parameters = {**parameters, 'output_dir': staging_dir}
        simulator(**staging_parameters)
        staged_raw = staging_dir / 'raw'
        if not staged_raw.is_dir():
            raise ValueError('simulator did not produce a raw directory')
        output_dir.mkdir(parents=True, exist_ok=True)
        for source in staged_raw.iterdir():
            shutil.move(str(source), str(output_dir / source.name))
        logger.success('Published dataset files to {}', output_dir)
