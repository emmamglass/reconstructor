#!/usr/bin/env
'''Reconstructor 

Reconstructor is an automatic genome scale metabolic network reconstruction tool that is user-friendly, COBRApy compatible, and uses a pFBA-based
gap-filling technique. 

Inputs
---------
Type 1: Annotated protein fasta file
Type 2: BLASTp output
Type 3: SBML Model

Output
---------
Well-annotated SBML model that uses ModelSEED namespace and is directly compatible with COBRApy without the need for additional compatibility modules

Example of how to run reconstructor
-----------------------------------
Type 1 input:  python -m reconstructor --input Osplanchnicus.aa.fasta --type 1 --gram negative --other_args <args>
Type 2 input: python -m reconstructor --input Osplanchnicus.hits.out --type 2 --gram negative --other_args <args>
Type 3 input: python -m reconstructor --input Osplanchnicus.sbml --type 3 --other_args <args>

Options for Running Reconstructor 
---------------------------------
--input <input file, Required>
--type <input file type, .fasta = 1, diamond blastp output = 2, .sbml = 3, Required, Default = 1> 
--gram <Type of Gram classificiation (positive or negative), default = positive>
--media <List of metabolites composing the media condition. Not required.>
--tasks <List of metabolic tasks. Not required>
--org <KEGG organism code. Not required>
--min_frac <Minimum objective fraction required during gapfilling, default = 0.01>
--max_frac <Maximum objective fraction allowed during gapfilling, default = 0.5>
--out <Name of output GENRE file, default = default>
--name <ID of output GENRE, default = default>
--cpu <Number of processors to use, default = 1>
--test <run installation tests, default = no>
'''

# Dependencies
import wget
import os
import cobra
import pickle
import argparse
import warnings
from multiprocessing import cpu_count
from sys import stdout
import platform
from unipath import Path
from cobra.manipulation.delete import *

from reconstructor._funcs import (
    OperatingSystemError,
    _run_blast,
    _read_blast,
    _genes_to_rxns,
    _create_model,
    _add_names,
    _find_reactions,
    _gapfill_model,
    _set_base_inputs,
    _add_annotation,
    _checkModel
)

# User defined arguments
parser = argparse.ArgumentParser(description='Generate genome-scale metabolic network reconstruction from KEGG BLAST hits.')
parser.add_argument('--input_file', default='none')
parser.add_argument('--file_type', default=1, help='Input file type: fasta=1, diamond blastp output=2, genre sbml=3')
parser.add_argument('--media', default='rich', help='List of metabolites composing the media condition. Not required.')
parser.add_argument('--tasks', default=[], help='List of metabolic tasks. Not required.')
parser.add_argument('--org', default='default', help='KEGG organism code. Not required.')
parser.add_argument('--min_frac', default=0.01, help='Minimum objective fraction required during gapfilling')
parser.add_argument('--max_frac', default=0.5, help='Maximum objective fraction allowed during gapfilling')
parser.add_argument('--gram', default='none', help='Type of Gram classificiation (positive or negative)')
parser.add_argument('--out', default='default', help='Name of output GENRE file')
parser.add_argument('--name', default='default', help='ID of output GENRE')
parser.add_argument('--cpu', default=1, help='Number of processors to use')
parser.add_argument('--gapfill', default='yes', help='gapfill your model?')
parser.add_argument('--exchange', default = 1, help='open exchange: 1, shut down exchange: 0')
parser.add_argument('--test', default = 'no', help='do you want to perform the test suite?')
args = parser.parse_args()



#----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = str(args.input_file)
    out_file = str(args.out)
    file_type = int(args.file_type)
    name = str(args.name)
    org = str(args.org)
    try:
        media = str(args.media).split(",")
        media = list(media)
    except:
        media = str(args.media)
    if media == ['default_media']:
        media = 'default'
    if media == ['rich']:
        media = 'rich'
    if media ==['minimal']:
        media = 'minimal'
    print(media)
    min_frac = float(args.min_frac)
    max_frac = float(args.max_frac)
    metabolic_tasks = list(args.tasks)
    new_id = str(args.name)
    gram_type = str(args.gram)
    processors = int(args.cpu)
    gapfill = str(args.gapfill)
    exchange_arg = int(args.exchange)
    test = str(args.test)

    #----------------------------------------------------------------------------------------------------------------------#
    if test == 'yes':
        script_path = str(os.path.dirname(os.path.realpath(__file__)))
        platform =  platform.system()
        print(platform)

        if platform == 'Darwin':
            if os.path.exists(script_path+'/refs') == False:
                    os.makedirs(script_path+'/refs')
            if os.path.exists(script_path+'/refs/gene_modelseed.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/gene_modelseed.pickle"
                wget.download(url, out = script_path+'/refs')
            if os.path.exists(script_path+'/refs/gene_names.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/gene_names.pickle"
                wget.download(url, out = script_path+'/refs')
            if os.path.exists(script_path+'/refs/screened_kegg_prokaryotes_pep_db.dmnd') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
                wget.download(url, out = script_path+'/refs')
            if os.path.exists(script_path+'/refs/universal.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/universal.pickle"
                wget.download(url, out = script_path+'/refs')

        if platform == 'Windows':
            if os.path.exists(script_path+r'\refs') == False:
                os.makedirs(script_path+r'\refs')
                print('script path', script_path)
            if os.path.exists(script_path+r'\refs\gene_modelseed.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/gene_modelseed.pickle"
                wget.download(url, out = script_path+r'\refs')
            if os.path.exists(script_path+r'\refs\gene_names.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/gene_names.pickle"
                wget.download(url, out = script_path+r'\refs')
            if os.path.exists(script_path+r'\refs\screened_kegg_prokaryotes_pep_db.dmnd') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
                wget.download(url, out = script_path+r'\refs')
            if os.path.exists(script_path+r'\refs\universal.pickle') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/universal.pickle"
                wget.download(url, out = script_path+r'\refs')

        if platform == 'Linux':
            raise OperatingSystemError("Reconstructor does not currently support Linux")

    #import test files 
    ##488.146.fa: an amino acid .fasta file used to test a type 1 input to reconstructor
    # JCP8151B.KEGGprot.out: a blast output file used to test a type 2 input to reconstructor
    # fmt.metaG.01044A.bin.149.KEGGprot.sbml: a .sbml genre used to test a type three input to reconstructor
        if platform == 'Darwin':
            if os.path.exists(script_path+'/testfiles') == False:
                os.makedirs(script_path+'/testfiles')
            if os.path.exists(script_path+'/testfiles/488.146.fa') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/488.146.fa"
                wget.download(url, out = script_path+'/testfiles')
            if os.path.exists(script_path+'/testfiles/JCP8151B.KEGGprot.out') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/JCP8151B.KEGGprot.out"
                wget.download(url, out = script_path+'/testfiles')
            if os.path.exists(script_path+'/testfiles/fmt.metaG.01044A.bin.149.KEGGprot.sbml') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/fmt.metaG.01044A.bin.149.KEGGprot.sbml"
                wget.download(url, out = script_path+'/testfiles')

        if platform == 'Windows':
            if os.path.exists(script_path+r'\testfiles' == False):
                os.makedirs(script_path+r'\testfiles')
            if os.path.exists(script_path+r'\testfiles\488.146.fa') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/488.146.fa"
                wget.download(url, out = script_path+r'\testfiles')
            if os.path.exists(script_path+r'\testfiles\JCP8151B.KEGGprot.out') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/JCP8151B.KEGGprot.out"
                wget.download(url, out = script_path+r'\testfiles')
            if os.path.exists(script_path+r'\testfiles\fmt.metaG.01044A.bin.149.KEGGprot.sbml') == False:
                url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/fmt.metaG.01044A.bin.149.KEGGprot.sbml"
                wget.download(url, out = script_path+r'\testfiles')

        exe = 'glpk_interface.py'
        '''if platform == 'Darwin' or platform == 'Linux':
            home_directory = os.path.expanduser('~')
        if platform == 'Windows':
            home_directory = os.path.expanduser(r'C:Users$USERNAME')'''
        if platform == 'Darwin':
            pa = Path(script_path).parent
            p = str(Path(pa).parent)
        if platform == 'Windows':
            p = str(Path(script_path).parent)
        #path = os.path.join(script_path, 'opt')
        for root, dirs, files in os.walk(p,topdown = True):
            for name in files:
                if name == exe:
                    file_path = os.path.abspath(os.path.join(root,name))
                    path = os.path.abspath(os.path.join(root))

                    if os.path.exists(file_path):
                        os.remove(file_path)
                
                    url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/glpk_interface.py"
                    wget.download(url, out = path)
        if platform == 'Darwin':
            cmd_line = "python -m reconstructor --input_file " + script_path+"/testfiles/488.146.fa --file_type 1 --gram negative"
            print("Performing test 1")
            print(cmd_line)
            os.system(cmd_line)

            cmd_line = "python -m reconstructor --input_file " + script_path+"/testfiles/JCP8151B.KEGGprot.out --file_type 2 --gram negative"
            print("Performing test 2")
            print(cmd_line)
            os.system(cmd_line)

            cmd_line = "python -m reconstructor --input_file " + script_path+"/testfiles/fmt.metaG.01044A.bin.149.KEGGprot.sbml --file_type 3 --gram negative"
            print("Performing test 3")
            print(cmd_line)
            os.system(cmd_line)

            quit()

        if platform == 'Windows':
            cmd_line = "python -m reconstructor --input_file " + script_path+r"\testfiles\488.146.fa --file_type 1 --gram negative"
            print("Performing test 1")
            print(cmd_line)
            os.system(cmd_line)

            cmd_line = "python -m reconstructor --input_file " + script_path+r"\testfiles\JCP8151B.KEGGprot.out --file_type 2 --gram negative"
            print("Performing test 2")
            print(cmd_line)
            os.system(cmd_line)

            cmd_line = "python -m reconstructor --input_file " + script_path+r"\testfiles\fmt.metaG.01044A.bin.149.KEGGprot.sbml --file_type 3 --gram negative"
            print("Performing test 3")
            print(cmd_line)
            os.system(cmd_line)

            quit()
    #----------------------------------------------------------------------------------------------------------------------#
    if gram_type == 'positive':
        print('\nUsing Gram positive objective function')
        universal_obj = 'biomass_GmPos'
    elif gram_type == 'negative':
        print('\nUsing Gram negative objective function')
        universal_obj = 'biomass_GmNeg'
    else:
        universal_obj = 'biomass'

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
    if media == 'rich':
        media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e','cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e'
    'cpd00322_e','cpd00066_e','cpd00054_e','cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e','cpd00104_e','cpd00149_e',
    'cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e','cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
    'cpd00092_e','cpd00117_e','cpd00067_e''cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e','cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e',
    'cpd00182_e','cpd00654_e','cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e','cpd00309_e','cpd00098_e','cpd00207_e',
    'cpd00082_e','cpd00129_e']
    elif media == 'minimal':
        media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e','cpd00084_e', 
    'cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e','cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e',
    'cpd00104_e','cpd00644_e','cpd00263_e','cpd00082_e']
    elif media == 'default':
        media = ['cpd00035_e','cpd00051_e','cpd00132_e','cpd00041_e','cpd00084_e','cpd00053_e','cpd00023_e',
    'cpd00033_e','cpd00119_e','cpd00322_e','cpd00107_e','cpd00039_e','cpd00060_e','cpd00066_e','cpd00129_e',
    'cpd00054_e','cpd00161_e','cpd00065_e','cpd00069_e','cpd00156_e','cpd00027_e','cpd00149_e','cpd00030_e',
    'cpd00254_e','cpd00971_e','cpd00063_e','cpd10515_e','cpd00205_e','cpd00099_e']
    else:
        media = media
    print(media)
    
    # Set media condition
    if len(media) != 0:
        media_condition = set(['EX_' + cpd for cpd in media])
        universal_reactions = set([x.id for x in universal.reactions])
        for rxn in universal_reactions:
            if rxn.startswith('EX_') == True:
                universal.reactions.get_by_id(rxn).bounds = (0, 1000.0)
            if rxn in media_condition:
                universal.reactions.get_by_id(rxn).bounds = (-1000.0, 10000)

    # Gapfill new model
    if gapfill == 'yes':
        if file_type != 3:
            print('Identifying new metabolism (Step 1 of 2)...')
        if file_type == 3:
            print('Identifying new metabolism...')
        draft_reactions = set([x.id for x in draft_genre.reactions])
        draft_metabolites = set([x.id for x in draft_genre.metabolites])
        warnings.filterwarnings('ignore')
        new_reactions = _find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
        print(new_reactions)
        filled_genre = _gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
        if file_type != 3:
            print('Identifying new metabolism (Step 2 of 2)...')
            filled_genre = _set_base_inputs(filled_genre, universal)
            media_reactions = _find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
            final_genre = _gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
            final_genre = _add_annotation(final_genre, gram_type)
        else: 
            final_genre = _add_annotation(filled_genre, universal_obj)
    else:
        draft_reactions = set([x.id for x in draft_genre.reactions])
        draft_metabolites = set([x.id for x in draft_genre.metabolites])
        final_genre = draft_genre
        final_genre = _add_annotation(final_genre, gram_type)
    warnings.filterwarnings('default')

    # Correct exchanges and check new model
    if exchange_arg == 0:
        for exch in final_genre.exchanges: exch.bounds = (0., 0.)
    else:
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

