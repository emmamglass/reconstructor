#!/usr/bin/env
'''USAGE: reconstruct.py input_file 
    Function that utilizes a pFBA solution with a universal reaction bag 
    in order to gapfill a model.
    
    Parameters
    ----------
    input_file : str
        Peptide fasta, diamond blastp output, or cobra model file
    type : int
        Input file type: fasta=1, diamond blastp output=2, genre sbml=3 (Default is 1)
    media : list
        List of strings of metabolite IDs that compase the media condition the new model should be gapfilled in (Default is complete media)
    tasks : list
        List of strings of reaction IDs that must be able to carry flux, not required (Default is none)
    min_frac : float
        Minimum objective flux required during gapfilling (Default is 0.1)
    max_frac : float
        Maximum objective flux allowed during gapfilling (Default is 1.0)
    gram : str
        Gram classificiation, positive or negative (Default is positive)
    out : str
        Name of output GENRE file in SBML format (Default is input_file + .sbml)
    id : str
        Name of output GENRE (Default is new_model)
    cpu : int
        Number of processors to use (Default is 1)
'''

# Dependencies
import os
import cobra
import pickle
import argparse
import warnings
import symengine
from random import shuffle
from multiprocessing import cpu_count
from sys import stdout
from copy import deepcopy
from subprocess import call
from cobra.util import solver
from cobra.manipulation.delete import *


# User defined arguments
parser = argparse.ArgumentParser(description='Generate genome-scale metabolic network reconstruction from KEGG BLAST hits.')
parser.add_argument('--input', default='none')
parser.add_argument('--type', default=1, help='Input file type: fasta=1, diamond blastp output=2, genre sbml=3')
parser.add_argument('--media', default=[], help='List of metabolites composing the media condition. Not required.')
parser.add_argument('--tasks', default=[], help='List of metabolic tasks. Not required.')
parser.add_argument('--org', default='default', help='KEGG organism code. Not required.')
parser.add_argument('--min_frac', default=0.01, help='Minimum objective fraction required during gapfilling')
parser.add_argument('--max_frac', default=0.5, help='Maximum objective fraction allowed during gapfilling')
parser.add_argument('--gram', default='positive', help='Type of Gram classificiation (positive or negative)')
parser.add_argument('--out', default='default', help='Name of output GENRE file')
parser.add_argument('--name', default='default', help='ID of output GENRE')
parser.add_argument('--cpu', default=1, help='Number of processors to use')
args = parser.parse_args()

# Run protein BLAST and save results
def _run_blast(inputfile, outputfile, database, processors, script_path):

    print('blasting %s vs %s'%(input_file,database))
    cmd_line = 'diamond blastp -p %s -d %s -q %s -o %s --more-sensitive --max-target-seqs 1 --quiet'%(processors,database,inputfile,outputfile)

    print('running blastp with the following command line...')
    print(cmd_line)
    os.system(cmd_line)
    print('finished blast')

    return outputfile

# Retreive KEGG hits
def _read_blast(blast_hits):
    hits = set()
    with open(blast_hits, 'r') as inFile:
        for line in inFile:
            line = line.split()
            hits |= set([line[1]])
    return hits

# Translate genes to ModelSEED reactions
def _genes_to_rxns(kegg_hits, gene_modelseed, organism):

    if organism != 'default':
        new_hits = _get_org_rxns(gene_modelseed, organism)
        gene_count = len(kegg_hits)
        kegg_hits |= new_hits
        gene_count = len(kegg_hits) - gene_count
        print('Added', gene_count, 'genes from', organism)

    rxn_db = {}
    for gene in kegg_hits:
        try:
            rxns = gene_modelseed[gene]
        except KeyError:
            continue

        for rxn in rxns:
            rxn = rxn + '_c'
            try:
                rxn_db[rxn].append(gene)
            except KeyError:
                rxn_db[rxn] = [gene]

    return rxn_db


# Get genes for organism from reference genome
def _get_org_rxns(gene_modelseed, organism):
    rxn_db = {}

    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(':')[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)


# Create draft GENRE and integrate GPRs
def _create_model(rxn_db, universal, input_id):
    new_model = cobra.Model('new_model')

    for x in rxn_db.keys():
        try:
            rxn = deepcopy(universal.reactions.get_by_id(x))
            new_model.add_reactions([rxn])
            new_model.reactions.get_by_id(x).gene_reaction_rule = ' or '.join(rxn_db[x])
        except KeyError:
            continue

    if input_id != 'default':
        new_model.id = input_id

    return new_model


# Add gene names
def _add_names(model, gene_db):
    for gene in model.genes:
        try:
            gene.name = gene_db[gene.id].title()
        except KeyError:
            continue

    return model


# pFBA gapfiller
def _find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type):

    stdout.write('\r[                                         ]')
    stdout.flush()

    # Modify universal reaction bag
    new_rxn_ids = set()
    with reaction_bag as universal:

        # Remove overlapping reactions from universal bag, and reset objective if needed
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()
        remove_rxns = []
        for rxn in model.reactions:
            if rxn.id == obj and file_type != 3:
                continue

            orig_rxn_ids |= set([rxn.id])
            try:
                test = universal.reactions.get_by_id(rxn.id)
                remove_rxns.append(rxn.id)
            except:
                continue

        # Add model reactions to universal bag
        universal.remove_reactions(list(set(remove_rxns)))
        add_rxns = []
        for x in model.reactions:
            if x.id != obj or file_type == 3:
                add_rxns.append(deepcopy(x))
        universal.add_reactions(add_rxns)

        # Set lower bounds for metaboloic tasks
        if len(tasks) != 0:
            for rxn in tasks:
                try:
                    universal.reactions.get_by_id(rxn).lower_bound = fraction
                except:
                    continue

        stdout.write('\r[---------------                          ]')
        stdout.flush()

        # Set minimum lower bound for previous objective
        universal.objective = universal.reactions.get_by_id(obj)           ############## problem here
        prev_obj_val = universal.slim_optimize()
        if step == 1:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
        	   lb=prev_obj_val*fraction, ub=prev_obj_val*max_fraction)
        elif step == 2:
            prev_obj_constraint = universal.problem.Constraint(universal.reactions.get_by_id(obj).flux_expression, 
               lb=prev_obj_val*max_fraction, ub=prev_obj_val)
        universal.solver.add(prev_obj_constraint)
        universal.solver.update()

        # Assemble forward and reverse components of all reactions
        coefficientDict = {}
        pfba_expr = symengine.RealDouble(0)
        for rxn in universal.reactions:
        	if rxn.id in orig_rxn_ids:
        		coefficientDict[rxn.forward_variable] = 0.0
        		coefficientDict[rxn.reverse_variable] = 0.0
        	else:
        		coefficientDict[rxn.forward_variable] = 1.0
        		coefficientDict[rxn.reverse_variable] = 1.0

        stdout.write('\r[--------------------------               ]')
        stdout.flush()

        # Create objective, based on pFBA
        universal.objective = 0
        universal.solver.update()
        universal.objective = universal.problem.Objective(symengine.RealDouble(0), direction='min', sloppy=True)
        universal.objective.set_linear_coefficients(coefficientDict)
        
        stdout.write('\r[----------------------------------       ]')
        stdout.flush()

        # Run FBA and identify reactions from universal that are now active
        solution = universal.optimize()
        
    new_rxn_ids = set([rxn.id for rxn in reaction_bag.reactions if abs(solution.fluxes[rxn.id]) > 1e-6]).difference(orig_rxn_ids)
    stdout.write('\r[-----------------------------------------]\n')
    warnings.filterwarnings('default')

    return(new_rxn_ids)    


# Add new reactions to model
def _gapfill_model(model, universal, new_rxn_ids, obj, step):

    # Get reactions and metabolites to be added to the model
    new_rxns = []
    if step == 1: new_rxns.append(deepcopy(universal.reactions.get_by_id(obj)))
    for rxn in new_rxn_ids: 
        if rxn != obj:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(rxn)))
    
    # Create gapfilled model 
    model.add_reactions(new_rxns)
    model.objective = model.problem.Objective(model.reactions.get_by_id(obj).flux_expression, direction='max')

    # Identify extracellular metabolites still need exchanges
    for cpd in model.metabolites:
        if cpd.compartment != 'extracellular':
            continue
        else:
            try:
                test = model.reactions.get_by_id('EX_' + cpd.id)
            except KeyError:
                exch_id = 'EX_' + cpd.id
                model.add_boundary(cpd, type='exchange', reaction_id=exch_id, lb=-1000.0, ub=1000.0)
                model.reactions.get_by_id(exch_id).name = cpd.name + ' exchange'

    return model


# Set uptake of specific metabolites in complete medium gap-filling
def _set_base_inputs(model, universal):
    tasks = ['EX_cpd00035_e','EX_cpd00051_e','EX_cpd00132_e','EX_cpd00041_e','EX_cpd00084_e','EX_cpd00053_e','EX_cpd00023_e',
    'EX_cpd00033_e','EX_cpd00119_e','EX_cpd00322_e','EX_cpd00107_e','EX_cpd00039_e','EX_cpd00060_e','EX_cpd00066_e','EX_cpd00129_e',
    'EX_cpd00054_e','EX_cpd00161_e','EX_cpd00065_e','EX_cpd00069_e','EX_cpd00156_e','EX_cpd00027_e','EX_cpd00149_e','EX_cpd00030_e',
    'EX_cpd00254_e','EX_cpd00971_e','EX_cpd00063_e','EX_cpd10515_e','EX_cpd00205_e','EX_cpd00099_e']

    new_rxns = []
    for exch in tasks: 
        try:
            test = model.reactions.get_by_id(exch)
        except:
            new_rxns.append(deepcopy(universal.reactions.get_by_id(exch)))
    model.add_reactions(new_rxns)
    for exch in tasks: model.reactions.get_by_id(exch).bounds = (-1000., -0.01)

    return model


def _add_annotation(model, obj='built'):
    
    # Genes
    for gene in model.genes:
        gene._annotation = {}
        gene.annotation['sbo'] = 'SBO:0000243'
        gene.annotation['kegg.genes'] = gene.id
    
    # Metabolites
    for cpd in model.metabolites: 
        cpd._annotation = {}
        cpd.annotation['sbo'] = 'SBO:0000247'
        if 'cpd' in cpd.id: cpd.annotation['seed.compound'] = cpd.id.split('_')[0]

    # Reactions
    for rxn in model.reactions:
        rxn._annotation = {}
        if 'rxn' in rxn.id: rxn.annotation['seed.reaction'] = rxn.id.split('_')[0]
        compartments = set([x.compartment for x in list(rxn.metabolites)])
        if len(list(rxn.metabolites)) == 1:
            rxn.annotation['sbo'] = 'SBO:0000627' # exchange
        elif len(compartments) > 1:
            rxn.annotation['sbo'] = 'SBO:0000185' # transport
        else:
            rxn.annotation['sbo'] = 'SBO:0000176' # metabolic

    # Biomass reactions
    if obj == 'built':
        model.reactions.EX_biomass.annotation['sbo'] = 'SBO:0000632'
        biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','peptidoglycan_rxn','lipid_rxn','cofactor_rxn','GmPos_cellwall','rxn10088_c','GmNeg_cellwall','biomass_rxn_gp','biomass_rxn_gn']
        for x in biomass_ids:
            try:
                model.reactions.get_by_id(x).annotation['sbo'] = 'SBO:0000629'
            except:
                continue
    else:
        model.reactions.get_by_id(obj).annotation['sbo'] = 'SBO:0000629'

    return model
    

# Run some basic checks on new models
def _checkModel(pre_reactions, pre_metabolites, post_model):

	# Check for objective flux
    new_genes = len(post_model.genes)
    new_rxn_ids = set([x.id for x in post_model.reactions]).difference(pre_reactions)
    new_cpd_ids = set([x.id for x in post_model.metabolites]).difference(pre_metabolites)
    test_flux = round(post_model.slim_optimize(), 3)

	
	# Report to user
    print('\tDraft reconstruction had', str(new_genes), 'genes,', str(len(pre_reactions)), 'reactions, and', str(len(pre_metabolites)), 'metabolites')
    print('\tGapfilled', str(len(new_rxn_ids)), 'reactions and', str(len(new_cpd_ids)), 'metabolites\n')
    print('\tFinal reconstruction has', str(len(post_model.reactions)), 'reactions and', str(len(post_model.metabolites)), 'metabolites')
    print('\tFinal objective flux is', str(round(test_flux, 3)))

    return test_flux

#----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = str(args.input)
    out_file = str(args.out)
    file_type = int(args.type)
    name = str(args.name)
    org = str(args.org)
    media = list(args.media)
    min_frac = float(args.min_frac)
    max_frac = float(args.max_frac)
    metabolic_tasks = list(args.tasks)
    new_id = str(args.name)
    gram_type = str(args.gram)
    processors = int(args.cpu)

    if gram_type == 'positive':
        print('\nUsing Gram positive objective function')
        universal_obj = 'biomass_GmPos'
    elif gram_type == 'negative':
        print('\nUsing Gram negative objective function')
        universal_obj = 'biomass_GmNeg'
    else:
        print('\nWARNING: Uknown Gram type selected. Using default of positive')
        universal_obj = 'biomass_GmPos'

    if min_frac <= 0.0 or min_frac > 1.0:
        print('WARNING: Improper minimum fraction selected. Defaulting to 0.01')
        min_frac = 0.01
    if max_frac <= 0.0 or max_frac > 1.0:
        print('WARNING: Improper maximum fraction selected. Defaulting to 0.5')
        max_frac = 0.5
    if max_frac < min_frac:
        print('WARNING: Input maximum fraction less than minimum fraction. Minimum set to half maximum')
        min_frac = max_frac * 0.5

    if org != 'default':
        print('Including additional genes from KEGG genome of', org)

    # Maximum fraction should not be too high, otherwise the gapfiller adds too many reactions
    print('Using minimum objective flux fraction of', min_frac,'and maximum fraction of', max_frac)

    if processors > cpu_count():
        print('WARNING: Requested more processors than are available. Using maximum of', cpu_count())
        processors = cpu_count()
    print('Using', processors, 'processor(s)\n')


    # Load databases
    print('Loading GENRE construction databases...')
    script_path = str(os.path.dirname(os.path.realpath(__file__)))
    kegg_prot_db = script_path + '/refs/screened_kegg_prokaryotes_pep_db'
    stdout.write('\r[                                         ]')
    stdout.flush()
    filename = script_path + '/refs/gene_modelseed.pickle'
    with open(filename, 'rb') as f: gene_modelseed = pickle.load(f)
    stdout.write('\r[---------------                          ]')
    stdout.flush()
    filename = script_path + '/refs/universal.pickle'
    with open(filename, 'rb') as f: universal = pickle.load(f)
    stdout.write('\r[------------------------------           ]')
    stdout.flush()
    filename = script_path + '/refs/gene_names.pickle'
    with open(filename, 'rb') as f: gene_names = pickle.load(f)
    stdout.write('\r[-----------------------------------------]\n')

    # Check input file type
    if file_type == 1:
        print('Aligning peptide sequences to KEGG database, may take some time...')
        blast_results = input_file.rstrip('fastn') + 'KEGGprot.out'
        print('Saving BLASTp results to', blast_results,'\n')
        _run_blast(input_file, blast_results, kegg_prot_db, str(processors), script_path)
    elif file_type == 2:
        blast_results = input_file
    else:
        try:
            draft_genre = cobra.io.read_sbml_model(input_file)
        except:
            draft_genre = cobra.io.load_json_model(input_file)


    # Handle gap-filling if that's all that is needed
    if file_type != 3:
        if blast_results != 'none':
            print('Creating draft GENRE from BLAST results...')
            gene_hits = _read_blast(blast_results)
        else: 
            gene_hits = set()
        rxns = _genes_to_rxns(gene_hits, gene_modelseed, org)
        draft_genre = _create_model(rxns, universal, new_id)
        draft_genre = _add_names(draft_genre, gene_names)
    else:
        universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]

    # Handle media conditions
    rich_media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e','cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e',
    'cpd00322_e','cpd00066_e','cpd00054_e','cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e','cpd00104_e','cpd00149_e',
    'cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e','cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
    'cpd00092_e','cpd00117_e','cpd00067_e''cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e','cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e',
    'cpd00182_e','cpd00654_e','cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e','cpd00309_e','cpd00098_e','cpd00207_e',
    'cpd00082_e','cpd00129_e']
    minimal_media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e','cpd00084_e', 
    'cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e','cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e',
    'cpd00104_e','cpd00644_e','cpd00263_e','cpd00082_e']
    if media == 'rich':
        media = rich_media
    elif media == 'minimal':
        media = minimal_media
    else:
        media = media
    # Set media condition
    if len(media) != 0:
        media_condition = set(['EX_' + cpd for cpd in media])
        universal_exchanges = set([x.id for x in universal.exchanges])
        for rxn in universal_exchanges:
            if rxn in media_condition:
                universal.reactions.get_by_id(rxn).bounds = (-1000.0, 1000.0)
            else:
                universal.reactions.get_by_id(rxn).bounds = (0.0, 1000.0)

    # Gapfill new model
    if file_type != 3:
        print('Identifying new metabolism (Step 1 of 2)...')
    if file_type == 3:
        print('Identifying new metabolism...')
    draft_reactions = set([x.id for x in draft_genre.reactions])
    draft_metabolites = set([x.id for x in draft_genre.metabolites])
    warnings.filterwarnings('ignore')
    new_reactions = _find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
    filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
    if file_type != 3:
        print('Identifying new metabolism (Step 2 of 2)...')
        filled_genre = _set_base_inputs(filled_genre, universal)
        media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
        final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
        final_genre = _add_annotation(final_genre)
    else: 
        final_genre = _add_annotation(filled_genre, universal_obj)
    warnings.filterwarnings('default')

    # Correct exchanges and check new model
    for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
    for rxn in final_genre.reactions:
        if 'Exchange reaction for' in rxn.name:
            rxn.name = list(rxn.metabolites)[0].name + ' exchange'
    biomass = _checkModel(draft_reactions, draft_metabolites, final_genre)

    # Write new model to sbml
    input_file = input_file.split('/')[-1] # write to working directory
    if file_type == 1:
        if new_id != 'default':
            out_file = input_file.rstrip('fastn') + new_id + '.sbml'
        else:
            out_file = input_file.rstrip('fastn') + 'sbml'
    elif file_type == 2:
        if new_id != 'default':
            if input_file != 'none':
                out_file = input_file.rstrip('out') + new_id + '.sbml'
            else:
                out_file = new_id + '.sbml'
        else:
            if org != 'default':
                out_file = org + '.sbml'
            else:
                out_file = input_file.rstrip('out') + 'sbml'
    elif file_type == 3:
        if new_id != 'default':
            out_file = input_file.rstrip('sbml') + new_id + '.extended.sbml'
        else:
            out_file = input_file.rstrip('sbml') + 'extended.sbml'

    print('\nSaving new GENRE to', out_file, '\n')
    cobra.io.write_sbml_model(final_genre, out_file)

