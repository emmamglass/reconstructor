from pathlib import Path
from tempfile import TemporaryDirectory as TempDir

import pytest

from reconstructor._funcs import run_blast, read_blast, genes_to_rxns


def test_blast(tiny_fasta_file: Path, expected_blast_output_file: Path, kegg_prot_db: Path):
    """
    
    """
    with TempDir() as tmpdir:
        out_path = Path(tmpdir).joinpath("blast.out")
        run_blast(tiny_fasta_file, out_path, kegg_prot_db, 1, None)
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
        "aai:AARI_04680",
        "gvg:HMPREF0421_20718",
        "cst:CLOST_0588",
        "vei:Veis_0353",
        "dwd:DSCW_47150",
        "gvh:HMPREF9231_1083",
        "jli:EXU32_08840",
        "gvh:HMPREF9231_0377",
        "cth:Cthe_0797",
        "bbp:BBPR_1794"
    }
    result = read_blast(blast_output_file)
    assert result == expected


def test_genes_to_rxns(blast_output_file: Path, modelseed_db: dict[str, list[str]]):
    """
    Simple test to ensure that a blast output file is correctly read and
    converted to a dictionary of reactions and associated genes.
    """
    expected = {
        'rxn38278_c': ['vei:Veis_0353'],
        'rxn32389_c': ['vei:Veis_0353'],
        'rxn38702_c': ['vei:Veis_0353'],
        'rxn03869_c': ['vei:Veis_0353'],
        'rxn06522_c': ['vei:Veis_0353'],
        'rxn38142_c': ['vei:Veis_0353'],
        'rxn02518_c': ['cst:CLOST_0588'],
        'rxn01962_c': ['cst:CLOST_0588'],
        'rxn33568_c': ['cst:CLOST_0588'],
        'rxn32974_c': ['cst:CLOST_0588'],
        'rxn32148_c': ['cst:CLOST_0588'],
        'rxn32257_c': ['cst:CLOST_0588'],
        'rxn33469_c': ['gvh:HMPREF9231_0377'],
        'rxn01987_c': ['gvh:HMPREF9231_0377'],
        'rxn02008_c': ['jli:EXU32_08840'],
        'rxn15597_c': ['cth:Cthe_0797'],
        'rxn30045_c': ['dwd:DSCW_47150'],
        'rxn16583_c': ['aai:AARI_04680'],
        'rxn15443_c': ['bbp:BBPR_1794'],
        'rxn37953_c': ['gvg:HMPREF0421_20718'],
        'rxn03408_c': ['gvh:HMPREF9231_1083'],
        'rxn03933_c': ['gvh:HMPREF9231_1083']
    }

    blast_hits = read_blast(blast_output_file)
    reactions = genes_to_rxns(blast_hits, modelseed_db, "default")
    assert reactions == expected
