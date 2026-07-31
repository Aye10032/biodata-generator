import tomllib
from dataclasses import dataclass
from pathlib import Path

PATH_PARAMETERS = frozenset(
    {
        'reference',
        'candidate_targets',
        'transcripts_path',
        'annotation_path',
        'output_dir',
    }
)


@dataclass(frozen=True)
class BatchCase:
    name: str
    simulator: str
    enabled: bool
    parameters: dict[str, object]


@dataclass(frozen=True)
class BatchConfig:
    source: Path
    workspace_root: Path
    cases: tuple[BatchCase, ...]


def load_batch_config(path: Path) -> BatchConfig:
    source = path.resolve()
    with source.open('rb') as handle:
        document = tomllib.load(handle)

    version = document.get('version')
    if version != 2:
        raise ValueError(f'{path}: unsupported or missing config version {version!r}')

    workspace_value = document.get('workspace_root', '.')
    if not isinstance(workspace_value, str):
        raise TypeError(f'{path}: workspace_root must be a string')
    workspace_root = _resolve_path(source.parent, workspace_value)

    defaults = document.get('defaults', {})
    if not isinstance(defaults, dict):
        raise TypeError(f'{path}: defaults must be a table')

    raw_cases = document.get('cases')
    if not isinstance(raw_cases, dict) or not raw_cases:
        raise ValueError(f'{path}: at least one [cases.<name>] entry is required')

    cases: list[BatchCase] = []
    for name, raw_case in raw_cases.items():
        if not isinstance(raw_case, dict):
            raise TypeError(f'{path}: case {name!r} must be a table')
        simulator = raw_case.get('simulator')
        enabled = raw_case.get('enabled', True)
        parameters = raw_case.get('parameters', {})
        if not isinstance(simulator, str) or not simulator:
            raise ValueError(f'{path}: case {name!r} has an invalid simulator')
        if not isinstance(enabled, bool):
            raise TypeError(f'{path}: case {name!r} enabled must be boolean')
        if not isinstance(parameters, dict):
            raise TypeError(f'{path}: case {name!r} parameters must be a table')

        merged: dict[str, object] = {**defaults, **parameters}
        for parameter in PATH_PARAMETERS:
            value = merged.get(parameter)
            if value is not None:
                if not isinstance(value, str):
                    raise ValueError(f'{path}: case {name!r} parameter {parameter!r} must be a path string')
                merged[parameter] = _resolve_path(workspace_root, value)
        if 'output_dir' not in merged:
            raise ValueError(f'{path}: case {name!r} is missing parameters.output_dir')

        cases.append(BatchCase(name, simulator, enabled, merged))

    return BatchConfig(source, workspace_root, tuple(cases))


def _resolve_path(base: Path, value: str) -> Path:
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return candidate.resolve()
