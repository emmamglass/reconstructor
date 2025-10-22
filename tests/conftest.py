from pathlib import Path
from importlib import resources

import pytest
import cobra

import reconstructor.resources


@pytest.fixture
def resource_dir() -> Path:
    return resources.files(__package__).joinpath("resources")


@pytest.fixture
def modelseed_db() -> dict[str, list[str]]:
    return reconstructor.resources.get_gene_mseed_map()
    

@pytest.fixture
def universal_model() -> cobra.Model:
    return reconstructor.resources.get_universal_model()


@pytest.fixture
def kegg_prot_db() -> Path:
    return reconstructor.resources.get_diamond_db_path()
    

@pytest.fixture
def blast_output_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("sample_blast.out")


@pytest.fixture
def tiny_fasta_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("tiny_fasta.fa")


@pytest.fixture
def expected_blast_output_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("tiny_blast_expected.out")
