from typing import Generator
from functools import lru_cache
from importlib import resources
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import contextmanager
import gzip
import json
import platform
import zipfile

import wget
import cobra


RESOURCE_DIR = resources.files(__package__)


@lru_cache(maxsize=None)
def get_universal_model() -> cobra.Model:
    resource = RESOURCE_DIR.joinpath("universal.sbml.gz")
    return cobra.io.read_sbml_model(resource)


@lru_cache(maxsize=None)
def get_gene_name_map() -> dict[str, str]:
    resource = RESOURCE_DIR.joinpath("gene_names.json.gz")
    with gzip.open(resource, "rt") as f:
        return json.load(f)


@lru_cache(maxsize=None)
def get_gene_mseed_map() -> dict[str, list[str]]:
    resource = RESOURCE_DIR.joinpath("gene_modelseed.json.gz")
    with gzip.open(resource, "rt") as f:
        return json.load(f)


def get_diamond_db_path() -> Path:
    return Path(RESOURCE_DIR.joinpath("screened_kegg_prokaryotes_pep_db.dmnd"))


def download_diamond_db():
    """
    Downloads the DIAMOND database file from the Reconstructor releases page.
    """
    url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
    wget.download(url, out=get_diamond_db_path())


@contextmanager
def diamond_exe() -> Generator[Path, None, None]:
    """
    Context manager returning a filepath to the DIAMOND executable.

    This involves extracting the appropriate executable from the zip archive
    containing all the DIAMOND executables into a temporary directory. The path
    to the executable in the resulting temporary directory is then returned.
    """
    with TemporaryDirectory() as tempdir:
        yield _unpack_diamond_exe(platform.system(), tempdir)


def _unpack_diamond_exe(system: str, dir = None) -> Path:
    """
    Unpacks the DIAMOND appropriate DIAMOND binary for the provided system into
    the specified directory and returns the path to the binary.
    """
    zip_path = RESOURCE_DIR.joinpath("diamond.zip")
    with zipfile.ZipFile(zip_path, "r", zipfile.ZIP_DEFLATED) as archive:
        pathstr = archive.extract(f"diamond-{system}", dir)
    return Path(pathstr)
