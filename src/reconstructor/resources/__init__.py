from functools import lru_cache
from importlib import resources
from pathlib import Path
import gzip
import json

import cobra

from reconstructor.utils import download, DownloadProgress


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


def download_diamond_db() -> Path:
    """
    Downloads the DIAMOND database file from the Reconstructor releases page.
    """
    url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
    progress = DownloadProgress("Downloading the DIAMOND database for blasting...")
    return download(url, path=get_diamond_db_path(), callback=progress)


def remove_diamond_db():
    """
    Deletes the DIAMOND database file. This is meant to be used prior to
    uninstalling Reconstructor; otherwise, the DIAMOND database will be left
    behind when Reconstructor is uninstalled.
    """
    path = get_diamond_db_path()
    path.unlink(missing_ok=True)
