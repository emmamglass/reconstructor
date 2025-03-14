import os
import platform
from copy import deepcopy
from sys import stdout
import warnings

import cobra
import symengine


class OperatingSystemError(Exception):
    pass


# Run protein BLAST and save results
def _run_blast(inputfile, outputfile, database, processors, script_path):
    ''' runs protein BLAST and saves results '''
    opsys =  platform.system()
    print(opsys)
    print('blasting %s vs %s'%(inputfile,database))
    if opsys == "Darwin":
        cmd_line = 'diamond blastp -p %s -d %s -q %s -o %s --more-sensitive --max-target-seqs 1 --quiet'%(processors,database,inputfile,outputfile)
        print('running blastp with the following command line...')
        print(cmd_line)
        os.system(cmd_line) 
    elif opsys == "Windows":
        script_path = str(os.path.dirname(os.path.realpath(__file__)))
        diamond_path = script_path + r"\diamond"
        cmd_line = str(diamond_path)+' blastp -p %s -d %s -q %s -o %s --more-sensitive --max-target-seqs 1 --quiet'%(processors,database,inputfile,outputfile)
        print('running blastp with the following command line...')
        print(cmd_line)
        os.system(cmd_line) 
    elif opsys == 'Linux':
        raise OperatingSystemError("Reconstructor does not currently support Linux")
    
    print('finished blast')

    return outputfile


# Retreive KEGG hits
def _read_blast(blast_hits):
    ''' Retrieves KEGG hits '''
    hits = set()
    with open(blast_hits, 'r') as inFile:
        for line in inFile:
            line = line.split()
            hits |= set([line[1]])
    return hits


# Translate genes to ModelSEED reactions
def _genes_to_rxns(kegg_hits, gene_modelseed, organism):
    ''' Translates genes to ModelSEED reactions '''
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
    ''' Get genes for organism from reference genome '''
    rxn_db = {}

    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(':')[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)


# Create draft GENRE and integrate GPRs
def _create_model(rxn_db, universal, input_id):
    ''' Create draft GENRE and integrate GPRs '''
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
    ''' Add gene names '''
    for gene in model.genes:
        try:
            gene.name = gene_db[gene.id].title()
        except KeyError:
            continue

    return model


# pFBA gapfiller
def _find_reactions(model, reaction_bag, tasks, obj, fraction, max_fraction, step, file_type):
    ''' pFBA gapfiller that modifies universal reaction bag, removes overlapping reacitons from universal reaction bag
    and resets objective if needed, adds model reaction to universal bag, sets lower bound for metabolic tasks, 
    sets minimum lower bound for previous objective, assemble forward and reverse components of all reactions,
    create objective, based on pFBA, run FBA and identify reactions from universal that are now active'''
    stdout.write('\r[                                         ]')
    stdout.flush()

    # Modify universal reaction bag
    new_rxn_ids = set() #make empty set we will add new reaction ids to
    with reaction_bag as universal: #set the reaction bag as the universal reaction databse

        # Remove overlapping reactions from universal bag, and reset objective if needed
        warnings.filterwarnings('ignore')
        orig_rxn_ids = set()    #original reaction ids start as an empty set
        remove_rxns = []    #reactions to remove is empty vector
        for rxn in model.reactions: #for a reaction in the draft model reactions
            if rxn.id == obj and file_type != 3: #if a reaction is part of the objective function 
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
        universal.objective = universal.reactions.get_by_id(obj) 
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
    '''Adds new reactions to model by getting reactions and metabolites to be added to the model, creates gapfilled model, 
    and identifies extracellular metabolites that still need exchanges '''
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


def _add_annotation(model, gram, obj='built'):
    ''' Add gene, metabolite, reaction ,biomass reaction annotations '''
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
        try:
            model.reactions.EX_biomass.annotation['sbo'] = 'SBO:0000632'
        except:
            pass
        if gram == 'none':
            biomass_ids = ['dna_rxn','rna_rxn','protein_rxn','teichoicacid_rxn','lipid_rxn','cofactor_rxn','rxn10088_c','biomass_rxn']
        else:
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
    ''' Run basic checks on new models (checking for objective flux'''

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
