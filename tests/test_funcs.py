from pathlib import Path
from tempfile import TemporaryDirectory as TempDir

import pytest

from reconstructor._funcs import run_blast, read_blast, genes_to_rxns


def test_blast(tiny_fasta_file: Path, expected_blast_output_file: Path, kegg_prot_db: Path):
    """
    
    """
    with TempDir() as tmpdir:
        out_path = Path(tmpdir).joinpath("blast.out")
        run_blast(tiny_fasta_file, out_path, kegg_prot_db, 1)
        with out_path.open("r") as f:
            results = f.read()
        with expected_blast_output_file.open("w") as f:
            f.write(results)


def test_read_blast(blast_output_file: Path):
    """
    This is just a very basic test to ensure that a BLAST output file is read
    and parsed correctly.
    """
    expected = {
        "WP_004111608.1": "aai:AARI_04680",
        "WP_004111823.1": "gvg:HMPREF0421_20718",
        "WP_004112121.1": "cst:CLOST_0588",
        "WP_004113321.1": "vei:Veis_0353",
        "WP_004113429.1": "dwd:DSCW_47150",
        "WP_196227846.1": "gvh:HMPREF9231_1083",
        "WP_196227848.1": "jli:EXU32_08840",
        "WP_196227849.1": "gvh:HMPREF9231_0377",
        "WP_230479112.1": "cth:Cthe_0797",
        "WP_230479416.1": "bbp:BBPR_1794"
    }
    result = read_blast(blast_output_file)
    assert result == expected


def test_genes_to_rxns(blast_output_file: Path, modelseed_db: dict[str, list[str]]):
    """
    Simple test to ensure that a blast output file is correctly read and
    converted to a dictionary of reactions and associated genes.
    """
    expected = {
        "rxn16583_c": ["WP_004111608.1"],
        "rxn37953_c": ["WP_004111823.1"],
        "rxn02518_c": ["WP_004112121.1"],
        "rxn01962_c": ["WP_004112121.1"],
        "rxn33568_c": ["WP_004112121.1"],
        "rxn32974_c": ["WP_004112121.1"],
        "rxn32148_c": ["WP_004112121.1"],
        "rxn32257_c": ["WP_004112121.1"],
        "rxn38278_c": ["WP_004113321.1"],
        "rxn32389_c": ["WP_004113321.1"],
        "rxn38702_c": ["WP_004113321.1"],
        "rxn03869_c": ["WP_004113321.1"],
        "rxn06522_c": ["WP_004113321.1"],
        "rxn38142_c": ["WP_004113321.1"],
        "rxn30045_c": ["WP_004113429.1"],
        "rxn03408_c": ["WP_196227846.1"],
        "rxn03933_c": ["WP_196227846.1"],
        "rxn02008_c": ["WP_196227848.1"],
        "rxn33469_c": ["WP_196227849.1"],
        "rxn01987_c": ["WP_196227849.1"],
        "rxn15597_c": ["WP_230479112.1"],
        "rxn15443_c": ["WP_230479416.1"]
    }

    blast_hits = read_blast(blast_output_file)
    reactions = genes_to_rxns(blast_hits, modelseed_db, "default")
    assert reactions == expected
