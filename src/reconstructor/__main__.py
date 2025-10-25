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
import os
from pathlib import Path
import argparse
from multiprocessing import cpu_count
from tempfile import TemporaryDirectory
import zipfile

import cobra

from reconstructor._funcs import (
    run_blast,
    read_blast,
    genes_to_rxns,
    create_model,
    add_names,
    find_reactions,
    gapfill_model,
    set_base_inputs,
    add_annotation,
    check_model
)
from reconstructor.diamond import Diamond, download_diamond, DEFAULT_DIAMOND_VERSION
from reconstructor import resources, errors


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

# Diamond download options (only used when running the test suite)
group = parser.add_mutually_exclusive_group(required=False)
group.add_argument(
    "--diamond",
    nargs="?",
    const=DEFAULT_DIAMOND_VERSION,
    default=None,
    help="Force DIAMOND to be downloaded when running the test suite, and optionally specify the version"
)
group.add_argument(
    "--skip-diamond",
    action="store_true",
    default=False,
    help="Skip downloading a DIAMOND binary if running the test suite"
)

args = parser.parse_args()



#----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = str(args.input_file)
    out_file = str(args.out)
    file_type = int(args.file_type)
    name = str(args.name)
    org = str(args.org)
    if args.media in ["default", "default_media", "rich", "minimal"]:
        media = str(args.media)
    else:
        media = str(args.media).split(",")
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

        # Download a diamond binary if needed
        if args.diamond is not None:
            print(f"Getting DIAMOND v{args.diamond} from https://github.com/bbuchfink/diamond/releases")
            download_diamond(diamond_version=args.diamond)
        if not args.skip_diamond:
            try:
                diamond = Diamond()
            except errors.DiamondNotFoundError:
                print("DIAMOND not found...getting DIAMOND from https://github.com/bbuchfink/diamond/releases")
                download_diamond()
                diamond = Diamond()
            finally:
                print(f"Using DIAMOND v{diamond.get_version()} at {diamond.path}")

        # Download the diamond database file if it hasn't been downloaded yet
        diamond_db_path = resources.get_diamond_db_path()
        if not diamond_db_path.exists():
            print("Downloading the DIAMOND database for blasting...")
            resources.download_diamond_db()
            print("Done")

        # Run the three tests (each with a different input file)
        # - 488.146.clean.fa: an amino acid .fasta file used to test a type 1 input to reconstructor
        # - JCP8151B.KEGGprot.out: a blast output file used to test a type 2 input to reconstructor
        # - fmt.metaG.01044A.bin.149.KEGGprot.sbml: a .sbml genre used to test a type three input to reconstructor
        test_file_names = ["488.146.clean.fa", "JCP8151B.KEGGprot.out", "fmt.metaG.01044A.bin.149.KEGGprot.sbml"]
        input_types = [1, 2, 3]
        test_num = 0
        for test_file_name, input_type in zip(test_file_names, input_types):
            test_num += 1
            print(f"Performing test {test_num}")
            
            # Temporary directory to hold test files (and clean them up after test finishes)
            with TemporaryDirectory(dir=resources.RESOURCE_DIR) as tempdir:

                # Extract the test file
                with zipfile.ZipFile(resources.RESOURCE_DIR.joinpath("testfiles.zip")) as archive:
                    test_file = archive.extract(test_file_name, tempdir)

                # Run reconstructor
                output_file = Path(test_file).with_suffix(".out.sbml")
                cmd = f"python -m reconstructor --input_file {test_file} --file_type {input_type} --out {output_file} --gram negative --cpu {processors}"
                print(cmd)
                exitcode = os.system(cmd)
                if exitcode != 0:
                    raise errors.ReconstructorError(f"Test {test_num} failed with an error")

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
    kegg_prot_db = resources.get_diamond_db_path()
    print('\r[                                         ]', end='', flush=True)
    gene_modelseed = resources.get_gene_mseed_map()
    print('\r[---------------                          ]', end='', flush=True)
    universal = resources.get_universal_model()
    print('\r[------------------------------           ]', end='', flush=True)
    gene_names = resources.get_gene_name_map()
    print('\r[-----------------------------------------]')

    # Check input file type
    if file_type == 1:
        print('Aligning peptide sequences to KEGG database, may take some time...')
        blast_results = input_file.rstrip('fastn') + 'KEGGprot.out'
        print('Blast results will be saved to', blast_results,'\n')
        run_blast(input_file, blast_results, kegg_prot_db, str(processors))
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
            gene_hits = read_blast(blast_results)
        else: 
            gene_hits = set()
        rxns = genes_to_rxns(gene_hits, gene_modelseed, org)
        draft_genre = create_model(rxns, universal, new_id)
        draft_genre = add_names(draft_genre, gene_names)
    else:
        universal_obj = str(draft_genre.objective.expression).split()[0].split('*')[-1]

    # Handle media conditions
    if media == 'rich':
        media = ['cpd00001_e','cpd00035_e','cpd00041_e','cpd00023_e','cpd00119_e','cpd00107_e','cpd00060_e','cpd00161_e','cpd00069_e','cpd00084_e','cpd00033_e',
    'cpd00322_e','cpd00066_e','cpd00054_e','cpd00065_e','cpd00156_e','cpd00220_e','cpd00644_e','cpd00393_e','cpd00133_e','cpd00263_e','cpd00104_e','cpd00149_e',
    'cpd00971_e','cpd00099_e','cpd00205_e','cpd00009_e','cpd00063_e','cpd00254_e','cpd10515_e','cpd00030_e','cpd00242_e','cpd00226_e','cpd01242_e','cpd00307_e',
    'cpd00092_e','cpd00117_e','cpd00067_e','cpd00567_e','cpd00132_e','cpd00210_e','cpd00320_e','cpd03279_e','cpd00246_e','cpd00311_e','cpd00367_e','cpd00277_e',
    'cpd00182_e','cpd00654_e','cpd00412_e','cpd00438_e','cpd00274_e','cpd00186_e','cpd00637_e','cpd00105_e','cpd00305_e','cpd00309_e','cpd00098_e','cpd00207_e',
    'cpd00082_e','cpd00129_e']
    elif media == 'minimal':
        media = ['cpd00001_e','cpd00065_e','cpd00060_e','cpd00322_e','cpd00129_e','cpd00156_e','cpd00107_e','cpd00084_e', 
    'cpd00149_e','cpd00099_e','cpd10515_e','cpd00030_e','cpd00254_e','cpd00063_e','cpd00205_e','cpd00009_e','cpd00971_e','cpd00242_e',
    'cpd00104_e','cpd00644_e','cpd00263_e','cpd00082_e']
    elif media == 'default' or media == "default_media":
        media = ['cpd00035_e','cpd00051_e','cpd00132_e','cpd00041_e','cpd00084_e','cpd00053_e','cpd00023_e',
    'cpd00033_e','cpd00119_e','cpd00322_e','cpd00107_e','cpd00039_e','cpd00060_e','cpd00066_e','cpd00129_e',
    'cpd00054_e','cpd00161_e','cpd00065_e','cpd00069_e','cpd00156_e','cpd00027_e','cpd00149_e','cpd00030_e',
    'cpd00254_e','cpd00971_e','cpd00063_e','cpd10515_e','cpd00205_e','cpd00099_e']
    else:
        media = media
    
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
        new_reactions = find_reactions(draft_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 1, file_type)
        print(new_reactions)
        filled_genre = gapfill_model(draft_genre, universal, new_reactions, universal_obj, 1)
        if file_type != 3:
            print('Identifying new metabolism (Step 2 of 2)...')
            filled_genre = set_base_inputs(filled_genre, universal)
            media_reactions = find_reactions(filled_genre, universal, metabolic_tasks, universal_obj, min_frac, max_frac, 2, file_type)
            final_genre = gapfill_model(filled_genre, universal, media_reactions, universal_obj, 2)
            final_genre = add_annotation(final_genre, gram_type)
        else: 
            final_genre = add_annotation(filled_genre, universal_obj)
    else:
        draft_reactions = set([x.id for x in draft_genre.reactions])
        draft_metabolites = set([x.id for x in draft_genre.metabolites])
        final_genre = draft_genre
        final_genre = add_annotation(final_genre, gram_type)

    # Correct exchanges and check new model
    if exchange_arg == 0:
        for exch in final_genre.exchanges: exch.bounds = (0., 0.)
    else:
        for exch in final_genre.exchanges: exch.bounds = (-1000., 1000.)
    for rxn in final_genre.reactions:
        if 'Exchange reaction for' in rxn.name:
            rxn.name = list(rxn.metabolites)[0].name + ' exchange'
    biomass = check_model(draft_reactions, draft_metabolites, final_genre)

    # Write new model to sbml
    if out_file == "default":
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

