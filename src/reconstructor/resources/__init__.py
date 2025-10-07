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
import stat

import wget
import cobra


RESOURCE_DIR = resources.files(__package__)


@lru_cache(maxsize=None)
def get_universal_model() -> cobra.Model:
    """
    Get the universal reaction model.

    The first time this function is called, the universal model is loaded from
    the resources directory (which can take 30+ seconds) and then is cached so
    that subsequent calls can simply return the model without loading it again.
    """
    resource = RESOURCE_DIR.joinpath("universal.sbml.gz")
    return cobra.io.read_sbml_model(resource)


@lru_cache(maxsize=None)
def get_gene_name_map() -> dict[str, str]:
    """
    Get the dictionary mapping KEGG gene IDs to gene names.
    """
    resource = RESOURCE_DIR.joinpath("gene_names.json.gz")
    with gzip.open(resource, "rt") as f:
        return json.load(f)


@lru_cache(maxsize=None)
def get_gene_mseed_map() -> dict[str, list[str]]:
    """
    Get the dictionary mapping KEGG gene IDs to ModelSEED reaction IDs.
    """
    resource = RESOURCE_DIR.joinpath("gene_modelseed.json.gz")
    with gzip.open(resource, "rt") as f:
        return json.load(f)


def get_diamond_db_path() -> Path:
    """
    Get the filepath to the KEGG peptide DIAMOND database for blasting.
    """
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
    with TemporaryDirectory(dir=RESOURCE_DIR) as tempdir:
        yield _unpack_diamond_exe(platform.system(), tempdir)


def _unpack_diamond_exe(system: str, dir = None) -> Path:
    """
    Unpacks the DIAMOND appropriate DIAMOND binary for the provided system into
    the specified directory and returns the path to the binary.
    """
    zip_path = RESOURCE_DIR.joinpath(f"diamond-{system}.zip")
    with zipfile.ZipFile(zip_path, "r", zipfile.ZIP_DEFLATED) as archive:
        pathstr = archive.extract(f"diamond-{system}", dir)
    path = Path(pathstr)

    # Make sure it is executable on MacOS or Linux
    if system in ["Darwin", "Linux"]:
        curr_mode = path.stat().st_mode
        new_mode = curr_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        path.chmod(new_mode)

    return path
