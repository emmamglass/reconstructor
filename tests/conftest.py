from pathlib import Path
from importlib import resources

import pytest
import cobra

import reconstructor.resources


@pytest.fixture
def base_dir() -> Path:
    return Path(__file__).parent.parent


@pytest.fixture
def resource_dir() -> Path:
    return resources.files(__package__).joinpath("resources")


@pytest.fixture
def modelseed_db() -> dict[str, list[str]]:
    reconstructor.resources.get_gene_mseed_map()
    

@pytest.fixture
def universal_model() -> cobra.Model:
    return reconstructor.resources.get_universal_model()


@pytest.fixture
def kegg_prot_db(base_dir: Path) -> Path:
    return base_dir / "src" / "reconstructor" / "refs" / "screened_kegg_prokaryotes_pep_db"
    

@pytest.fixture
def blast_output_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("sample_blast.out")


@pytest.fixture
def tiny_fasta_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("tiny_fasta.fa")


@pytest.fixture
def expected_blast_output_file(resource_dir: Path) -> Path:
    return resource_dir.joinpath("tiny_blast_expected.out")
