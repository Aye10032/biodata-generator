from pathlib import Path

from click.testing import CliRunner

from bioflow_sim.cli import cli
from bioflow_sim.core.config import load_batch_config


def write_batch_fixture(tmp_path: Path) -> Path:
    (tmp_path / 'reference.fa').write_text('>chr1\n' + 'ACGT' * 500 + '\n', encoding='utf-8')
    config = tmp_path / 'batch.toml'
    config.write_text(
        """
version = 2
workspace_root = "."

[defaults]
seed = 7

[cases.small_dna]
simulator = "dna"
enabled = true

[cases.small_dna.parameters]
reference = "reference.fa"
output_dir = "output/small_dna"
sample = "BATCH"
technology_name = "illumina-pe"
reads = 5
read_length = 50
fragment_mean = 150
fragment_sd = 10
long_read_mean = 500
long_read_sd = 20
snvs = 2

[cases.disabled_dna]
simulator = "dna"
enabled = false

[cases.disabled_dna.parameters]
reference = "reference.fa"
output_dir = "output/disabled_dna"
sample = "DISABLED"
technology_name = "illumina-se"
reads = 1
read_length = 50
fragment_mean = 150
fragment_sd = 10
long_read_mean = 500
long_read_sd = 20
snvs = 0
""".strip()
        + '\n',
        encoding='utf-8',
    )
    return config


def test_batch_list_dry_run_and_generation(tmp_path: Path) -> None:
    config = write_batch_fixture(tmp_path)
    runner = CliRunner()

    listed = runner.invoke(cli, ['batch', str(config), '--list-cases'])
    assert listed.exit_code == 0, listed.output
    assert 'small_dna\tdna\tenabled' in listed.output
    assert 'disabled_dna\tdna\tdisabled' in listed.output

    dry_run = runner.invoke(cli, ['batch', str(config), '--dry-run'])
    assert dry_run.exit_code == 0, dry_run.output
    assert 'planned=1' in dry_run.output
    assert not (tmp_path / 'output').exists()

    generated = runner.invoke(cli, ['batch', str(config), '--case', 'small_dna'])
    assert generated.exit_code == 0, generated.output
    assert 'completed=1' in generated.output
    assert (tmp_path / 'output/small_dna/raw/BATCH_R1.fastq.gz').exists()
    assert not (tmp_path / 'output/disabled_dna').exists()


def test_default_config_resolves_repository_references() -> None:
    config_paths = (
        Path('config/lambda.toml'),
        Path('config/ecoli.toml'),
        Path('config/yeast.toml'),
        Path('config/homo_chr21.toml'),
    )
    for config_path in config_paths:
        config = load_batch_config(config_path)
        assert config.cases
        for case in config.cases:
            for parameter in ('reference', 'candidate_targets', 'transcripts_path', 'annotation_path'):
                path = case.parameters.get(parameter)
                if path is not None:
                    assert isinstance(path, Path)
                    assert path.exists()
