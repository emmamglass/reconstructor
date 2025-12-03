from typing import Optional, Union, Sequence, Literal
from collections import defaultdict
import os

import cobra
from optlang.symbolics import Zero

from reconstructor.diamond import Diamond
from reconstructor import medium, utils
from reconstructor._version import __version__


def run_blast(
    inputfile: Union[str, os.PathLike],
    outputfile: Union[str, os.PathLike],
    database: Union[str, os.PathLike],
    processors: Optional[int] = None,
) -> Union[str, os.PathLike]:
    """
    Runs protein BLAST and saves results.
    """

    if processors is None:
        options = ["--more-sensitive", "--max-target-seqs", "1"]
    else:
        options = ["-p", processors, "--more-sensitive", "--max-target-seqs", "1"]

    diamond = Diamond()
    diamond.blastp(database, inputfile, outputfile, *options, capture_output=True)

    return outputfile


def read_blast(blast_hits: Union[str, os.PathLike]) -> list[cobra.Gene]:
    """
    Retrieves KEGG hits from blast output.
    """

    hits: list[cobra.Gene] = []
    with open(blast_hits, "r") as file:
        for line in file:
            query_id, kegg_id, *_ = line.split()
            gene = cobra.Gene(utils.sanitize_sbml_id(query_id))
            gene.annotation["kegg.genes"] = kegg_id
            hits.append(gene)
    return hits


def genes_to_rxns(
    kegg_hits: list[cobra.Gene],
    gene_modelseed: dict[str, list[str]],
    organism: Optional[str],
) -> defaultdict[str, list[cobra.Gene]]:
    """
    Translates genes to ModelSEED reactions
    """
    org_genes = set()
    if organism is not None:
        blasted_genes = set(g.annotation["kegg.genes"] for g in kegg_hits)
        org_genes = _get_org_rxns(gene_modelseed, organism).difference(blasted_genes)
        print(f"Adding {len(org_genes)} from organism {organism}")

    rxn_db: defaultdict[str, list[cobra.Gene]] = defaultdict(list)
    for gene in kegg_hits:
        for rxn in gene_modelseed.get(gene.annotation["kegg.genes"], []):
            rxn = rxn + "_c"
            rxn_db[rxn].append(gene)

    for kegg_gene in org_genes:
        for rxn in gene_modelseed.get(kegg_gene, []):
            rxn = rxn + "_c"
            gene = cobra.Gene(kegg_gene)
            gene.annotation["kegg.genes"] = kegg_gene
            rxn_db[rxn].append(gene)

    return rxn_db


def _get_org_rxns(gene_modelseed: dict[str, list[str]], organism: str) -> set[str]:
    """
    Get genes for organism from reference genome.
    """

    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(":")[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)


def create_model(
    rxn_db: dict[str, list[cobra.Gene]],
    universal: cobra.Model,
    model_id: Optional[str] = None,
) -> cobra.Model:
    """
    Create draft GENRE and integrate GPRs.
    """
    new_model = cobra.Model(model_id)
    new_model.notes["source"] = f"Reconstructor v{__version__}"

    # Add the genes to the model
    for gene_lst in rxn_db.values():
        for gene in gene_lst:
            if not new_model.genes.has_id(gene.id):
                new_model.genes.add(gene)

    # Add the reactions to the model
    for rxn_id in rxn_db.keys():
        if universal.reactions.has_id(rxn_id):
            rxn: cobra.Reaction = universal.reactions.get_by_id(rxn_id)
            new_model.add_reactions([rxn.copy()])
            gpr = " or ".join(g.id for g in rxn_db[rxn_id])
            new_model.reactions.get_by_id(rxn_id).gene_reaction_rule = gpr

    return new_model


def add_names(model: cobra.Model, gene_db: dict[str, str]) -> cobra.Model:
    """
    Add gene names.
    """

    # FIXME: the current gene database doesn't have any actual gene names in it
    for gene in model.genes:
        if gene.id in gene_db:
            gene.name = gene_db[gene.id].title()

    return model


def find_reactions(
    model: cobra.Model,
    reaction_bag: cobra.Model,
    tasks: Optional[Sequence[str]],
    obj: str,
    fraction: float,
    max_fraction: float,
    step: Literal[1, 2],
    file_type: Literal[1, 2, 3],
) -> set[str]:
    """
    pFBA gapfiller.

    Modifies universal reaction bag, removes overlapping reacitons from
    universal reaction bag and resets objective if needed, adds model reaction
    to universal bag, sets lower bound for metabolic tasks, sets minimum lower
    bound for previous objective, assemble forward and reverse components of all
    reactions, create objective, based on pFBA, run FBA and identify reactions
    from universal that are now active.
    """

    print("\r[                                         ]", end="", flush=True)

    # Modify universal reaction bag
    new_rxn_ids = set()  # make empty set we will add new reaction ids to
    with reaction_bag as universal:  # set the reaction bag as the universal reaction database

        # Remove overlapping reactions from universal bag, and reset objective if needed
        orig_rxn_ids = set()
        remove_rxns = set()
        for rxn in model.reactions:

            # skip reaction if it is the objective function
            if rxn.id == obj and file_type != 3:
                continue

            orig_rxn_ids.add(rxn.id)
            if universal.reactions.has_id(rxn.id):
                remove_rxns.add(rxn.id)

        # Add model reactions to universal bag
        universal.remove_reactions(remove_rxns)
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(x.copy())
        universal.add_reactions(add_rxns)

        # Set lower bounds for metabolic tasks
        if tasks is not None:
            for rxn in tasks:
                if universal.reactions.has_id(rxn):
                    universal.reactions.get_by_id(rxn).lower_bound = fraction

        print("\r[---------------                          ]", end="", flush=True)

        # Set minimum lower bound for previous objective
        universal.objective = universal.reactions.get_by_id(obj)
        prev_obj_val = universal.slim_optimize()
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(
                universal.reactions.get_by_id(obj).flux_expression,
                lb=prev_obj_val * fraction,
                ub=prev_obj_val * max_fraction,
            )
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(
                universal.reactions.get_by_id(obj).flux_expression,
                lb=prev_obj_val * max_fraction,
                ub=prev_obj_val,
            )
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        coefficientDict = {}
        for rxn in universal.reactions:
            if rxn.id in orig_rxn_ids:
                coefficientDict[rxn.forward_variable] = 0.0
                coefficientDict[rxn.reverse_variable] = 0.0
            else:
                coefficientDict[rxn.forward_variable] = 1.0
                coefficientDict[rxn.reverse_variable] = 1.0

        print("\r[--------------------------               ]", end="", flush=True)

        # Create objective, based on pFBA
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(
            Zero, direction="min", sloppy=True
        )
        universal.objective.set_linear_coefficients(coefficientDict)

        print("\r[----------------------------------       ]", end="", flush=True)

        # Run FBA and identify reactions from universal that are now active
        solution = universal.optimize()

    active_rxns = solution.fluxes[solution.fluxes.abs() > 1e-6]
    new_rxn_ids = set(active_rxns.index).difference(orig_rxn_ids)
    print("\r[-----------------------------------------]")

    return new_rxn_ids


def gapfill_model(
    model: cobra.Model,
    universal: cobra.Model,
    new_rxn_ids: Sequence[str],
    obj: str,
    step: Literal[1, 2],
) -> cobra.Model:
    """
    Adds new reactions to model by getting reactions and metabolites to be added
    to the model, creates gapfilled model, and identifies extracellular
    metabolites that still need exchanges.
    """

    # Get reactions and metabolites to be added to the model
    new_rxns = []
    if step == 1:
        new_rxns.append(universal.reactions.get_by_id(obj).copy())
    for rxn in new_rxn_ids:
        if rxn != obj:
            new_rxns.append(universal.reactions.get_by_id(rxn).copy())

    # Create gapfilled model
    model.add_reactions(new_rxns)
    model.objective = model.problem.Objective(
        model.reactions.get_by_id(obj).flux_expression, direction="max"
    )

    # Identify extracellular metabolites still need exchanges
    for cpd in model.metabolites:
        if cpd.compartment != "extracellular":
            continue
        else:
            exch_id = "EX_" + cpd.id
            if not model.reactions.has_id(exch_id):
                model.add_boundary(
                    cpd, type="exchange", reaction_id=exch_id, lb=-1000.0, ub=1000.0
                )
                model.reactions.get_by_id(exch_id).name = cpd.name + " exchange"

    return model


def set_base_inputs(model: cobra.Model, universal: cobra.Model) -> cobra.Model:
    """
    Set uptake of specific metabolites in complete medium gap-filling.
    """

    tasks = [f"EX_{cpd}" for cpd in medium.COMPLETE]

    new_rxns = []
    for exch in tasks:
        if not model.reactions.has_id(exch):
            new_rxns.append(universal.reactions.get_by_id(exch).copy())
    model.add_reactions(new_rxns)
    for exch in tasks:
        model.reactions.get_by_id(exch).bounds = (-1000.0, -0.01)

    return model


def add_annotation(model: cobra.Model, gram: str, obj="built") -> cobra.Model:
    """
    Add gene, metabolite, reaction, and biomass reaction annotations.
    """

    # Genes
    for gene in model.genes:
        gene.annotation["sbo"] = "SBO:0000243"

    # Metabolites
    for cpd in model.metabolites:
        cpd.annotation["sbo"] = "SBO:0000247"
        if "cpd" in cpd.id:
            cpd.annotation["seed.compound"] = cpd.id.split("_")[0]

    # Reactions
    for rxn in model.reactions:
        if "rxn" in rxn.id:
            rxn.annotation["seed.reaction"] = rxn.id.split("_")[0]
        compartments = set([x.compartment for x in list(rxn.metabolites)])
        if len(list(rxn.metabolites)) == 1:
            rxn.annotation["sbo"] = "SBO:0000627"  # exchange
        elif len(compartments) > 1:
            rxn.annotation["sbo"] = "SBO:0000185"  # transport
        else:
            rxn.annotation["sbo"] = "SBO:0000176"  # metabolic

    # Biomass reactions
    if obj == "built":
        try:
            model.reactions.EX_biomass.annotation["sbo"] = "SBO:0000632"
        except AttributeError:
            pass
        if gram == "none":
            biomass_ids = [
                "dna_rxn",
                "rna_rxn",
                "protein_rxn",
                "teichoicacid_rxn",
                "lipid_rxn",
                "cofactor_rxn",
                "rxn10088_c",
                "biomass_rxn",
            ]
        else:
            biomass_ids = [
                "dna_rxn",
                "rna_rxn",
                "protein_rxn",
                "teichoicacid_rxn",
                "peptidoglycan_rxn",
                "lipid_rxn",
                "cofactor_rxn",
                "GmPos_cellwall",
                "rxn10088_c",
                "GmNeg_cellwall",
                "biomass_rxn_gp",
                "biomass_rxn_gn",
            ]
        for x in biomass_ids:
            try:
                model.reactions.get_by_id(x).annotation["sbo"] = "SBO:0000629"
            except KeyError:
                continue
    else:
        model.reactions.get_by_id(obj).annotation["sbo"] = "SBO:0000629"

    return model


def check_model(
    pre_reactions: Sequence[str],
    pre_metabolites: Sequence[str],
    post_model: cobra.Model,
):
    """
    Run basic checks on new models (checking for objective flux).
    """

    # Check for objective flux
    new_genes = len(post_model.genes)
    new_rxn_ids = set([x.id for x in post_model.reactions]).difference(pre_reactions)
    new_cpd_ids = set([x.id for x in post_model.metabolites]).difference(
        pre_metabolites
    )
    test_flux = round(post_model.slim_optimize(), 3)

    # Report to user
    print(
        f"\tDraft reconstruction had {new_genes} genes, {len(pre_reactions)} reactions, and {len(pre_metabolites)} metabolites"
    )
    print(
        f"\tGapfilled {len(new_rxn_ids)} reactions and {len(new_cpd_ids)} metabolites\n"
    )
    print(
        f"\tFinal reconstruction has {len(post_model.reactions)} reactions and {len(post_model.metabolites)} metabolites"
    )
    print(f"\tFinal objective flux is {round(test_flux, 3)}")

    return test_flux
