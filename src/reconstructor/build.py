from typing import Union, Literal, Sequence, Optional
import warnings
import os
import enum
from pathlib import Path

import cobra
from cobra.io.sbml import CobraSBMLError

from reconstructor import resources, _funcs
from reconstructor.medium import get_medium


class FileType(enum.Enum):
    FASTA = 1
    BLAST_OUTPUT = 2
    SBML = 3


def reconstruct(
    input_file: Union[str, os.PathLike],
    file_type: Literal[1, 2, 3] = 1,
    media: Union[Literal["rich", "minimal", "complete"], Sequence[str]] = "rich",
    tasks: Optional[Sequence[str]] = None,
    org: Optional[str] = None,
    min_frac: float = 0.01,
    max_frac: float = 0.5,
    gram: Literal["positive", "negative"] = "positive",
    model_id: Optional[str] = None,
    cpu: Optional[int] = None,
    gapfill: bool = True,
    open_exchanges: bool = True,
) -> cobra.Model:
    """
    Reconstruct a genome-scale metabolic model from an input file.

    For a FASTA input, uses DIAMOND to find similar protein sequences in KEGG.
    The blast output is then used to create a draft metabolic model, which is
    gap-filled with a pFBA-based approach.

    Parameters
    ----------
    input_file : str | PathLike
        The file path to the input file. The input file should be 1) a FASTA
        with amino acid sequences, 2) a tsv output from DIAMOND blast, or 3) an
        SBML model file.

    file_type : Literal[1, 2, 3]
        1 for a FASTA file, 2 for a DIAMOND blast output file, or 3 for an SBML
        model file. The default is 1.

    media : str | Sequence[str]
        Either the name of a registered media composition ("rich", "minimal", or
        "complete") or a sequence of extracellular metabolite IDs. The default
        is "rich".

    tasks : Sequence[str], optional
        A sequence of reaction IDs that should be active when gap-filling.
        Default is None.

    org : str, optional
        A KEGG organism ID. If supplied, all genes from this organism will added
        to the DIAMOND blast results. Default is None.

    min_frac : float
        The minimum objective fraction to use for pFBA gap-filling. Default is
        0.01.

    max_frac : float
        The maximum objective fraction to use for pFBA gap-filling. Default is
        0.5.

    gram : str
        The gram stain of the bacteria for which a model is being reconstructed.
        This will be used to determine the objective function for the model.
        Should either be "positive" or "negative". Default is "positive".

    model_id : str, optional
        The ID to use for the model that is created.

    cpu : int, optional
        How many processors to use for DIAMOND. Default is None, which means
        that DIAMOND will automatically detect how many to use.

    gapfill : bool
        Whether or not to perform pFBA gap-filling on the model. Default is
        True.

    open_exchanges : bool
        Whether or not to set the exchanges of the returned model to open or
        not. Default is True.

    Returns
    -------
    cobra.Model
        The reconstructed metabolic model as a COBRApy model.
    """

    # Convert input_file to path and ensure it exists
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(
            f"The provided input file path does not exist: {input_file}"
        )

    # Convert file type to enum value
    ftype = FileType(file_type)

    # Select the biomass objective function
    if gram == "positive":
        universal_obj = "biomass_GmPos"
    elif gram == "negative":
        universal_obj = "biomass_GmNeg"
    else:
        raise ValueError(
            f"Unrecognized gram stain value: {repr(gram)}. Should be 'positive', 'negative', or None."
        )

    # Ensure min/max fractions are properly set
    if min_frac <= 0.0 or min_frac > 1.0:
        warnings.warn(
            f"Improper minimum fraction selected ({repr(min_frac)}). Defaulting to 0.01"
        )
        min_frac = 0.01
    if max_frac <= 0.0 or max_frac > 1.0:
        warnings.warn(
            f"Improper maximum fraction selected ({repr(max_frac)}). Defaulting to 0.5"
        )
        max_frac = 0.5
    if max_frac < min_frac:
        warnings.warn(
            f"Input maximum fraction ({repr(max_frac)}) less than minimum fraction ({repr(min_frac)}). Minimum set to half maximum"
        )
        min_frac = max_frac * 0.5

    # Determine the media for gap-filling
    if isinstance(media, str):
        media = get_medium(media)

    # Load universal reaction model
    print("Loading universal model...")
    universal = resources.get_universal_model()

    if ftype in [FileType.FASTA, FileType.BLAST_OUTPUT]:
        # Run BLAST if needed
        if ftype == FileType.FASTA:
            print(
                "Blasting input sequences against KEGG database, may take some time..."
            )
            blast_db = resources.get_diamond_db_path()
            blast_file = input_path.with_suffix(".KEGGprot.out")
            _funcs.run_blast(input_path, blast_file, blast_db, cpu)
        else:
            blast_file = input_path

        # Draft model from BLAST output
        print("Drafting model from BLAST results...")
        gene_hits = _funcs.read_blast(blast_file)

        gene_modelseed = resources.get_gene_mseed_map()
        rxns = _funcs.genes_to_rxns(gene_hits, gene_modelseed, org)

        draft_model = _funcs.create_model(rxns, universal, model_id)

        gene_names = resources.get_gene_name_map()
        draft_model = _funcs.add_names(draft_model, gene_names)

    else:
        try:
            draft_model = cobra.io.read_sbml_model(input_path)
        except CobraSBMLError:
            draft_model = cobra.io.load_json_model(input_path)
        universal_obj = str(draft_model.objective.expression).split()[0].split("*")[-1]

    # Reactions/metabolites in draft model
    draft_reactions = set([x.id for x in draft_model.reactions])
    draft_metabolites = set([x.id for x in draft_model.metabolites])

    if gapfill:
        print("Gapfilling model...")

        # Set the medium on the universal model
        medium_dict = {}
        for cpd in media:
            exch_id = f"EX_{cpd}"
            if universal.reactions.has_id(exch_id):
                medium_dict[exch_id] = 1000
        universal.medium = medium_dict

        # First round of gap-filling
        new_reactions = _funcs.find_reactions(
            draft_model,
            universal,
            tasks,
            universal_obj,
            min_frac,
            max_frac,
            1,
            file_type,
        )
        filled_model = _funcs.gapfill_model(
            draft_model, universal, new_reactions, universal_obj, 1
        )

        if ftype != FileType.SBML:
            filled_model = _funcs.set_base_inputs(filled_model, universal)
            media_reactions = _funcs.find_reactions(
                filled_model,
                universal,
                tasks,
                universal_obj,
                min_frac,
                max_frac,
                2,
                file_type,
            )
            final_model = _funcs.gapfill_model(
                filled_model, universal, media_reactions, universal_obj, 2
            )
            final_model = _funcs.add_annotation(final_model, gram)
        else:
            final_model = _funcs.add_annotation(filled_model, universal_obj)

    else:
        final_model = _funcs.add_annotation(draft_model, gram)

    exchange_bounds = (-1000.0, 1000.0) if open_exchanges else (0.0, 0.0)
    for exch in final_model.exchanges:
        exch.bounds = exchange_bounds

    # Check the final model
    _funcs.check_model(draft_reactions, draft_metabolites, final_model)

    return final_model
