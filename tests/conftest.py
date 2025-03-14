from pathlib import Path
import pickle
from importlib.resources import files

import pytest
import cobra


@pytest.fixture
def base_dir() -> Path:
    return Path(__file__).parent.parent


@pytest.fixture
def resource_dir() -> Path:
    return files(__package__).joinpath("resources")


@pytest.fixture
def modelseed_db(base_dir: Path) -> dict[str, list[str]]:
    db_path = base_dir / "src" / "reconstructor" / "refs" / "gene_modelseed.pickle"
    with db_path.open("rb") as db:
        return pickle.load(db)
    

@pytest.fixture
def universal_model(base_dir: Path) -> cobra.Model:
    model_path = base_dir / "src" / "reconstructor" / "refs" / "universal.pickle"
    with model_path.open("rb") as db:
        return pickle.load(db)


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
