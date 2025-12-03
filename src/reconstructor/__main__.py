"""
Reconstructor
-------------

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
"""

# Dependencies
from pathlib import Path
import argparse
from tempfile import TemporaryDirectory
import zipfile

import cobra

from reconstructor.diamond import (
    Diamond,
    download_diamond,
    DEFAULT_DIAMOND_VERSION,
    get_diamond_path,
)
from reconstructor import resources, reconstruct


# User defined arguments
parser = argparse.ArgumentParser(
    description="Generate genome-scale metabolic network reconstruction from KEGG BLAST hits."
)
parser.add_argument("--input_file", default=None)
parser.add_argument(
    "--file_type",
    default=1,
    help="Input file type: fasta=1, diamond blastp output=2, genre sbml=3",
)
parser.add_argument(
    "--media",
    default="rich",
    help="List of metabolites composing the media condition. Not required.",
)
parser.add_argument(
    "--tasks", default=None, help="List of metabolic tasks. Not required."
)
parser.add_argument("--org", default=None, help="KEGG organism code. Not required.")
parser.add_argument(
    "--min_frac",
    default=0.01,
    help="Minimum objective fraction required during gapfilling",
)
parser.add_argument(
    "--max_frac",
    default=0.5,
    help="Maximum objective fraction allowed during gapfilling",
)
parser.add_argument(
    "--gram",
    default="positive",
    help="Type of Gram classificiation (positive or negative)",
)
parser.add_argument("--out", default=None, help="Name of output GENRE file")
parser.add_argument("--name", default=None, help="ID of output GENRE")
parser.add_argument("--cpu", default=None, help="Number of processors to use")
parser.add_argument("--gapfill", default="yes", help="gapfill your model?")
parser.add_argument(
    "--exchange", default=1, help="open exchange: 1, shut down exchange: 0"
)
parser.add_argument(
    "--test", default="no", help="do you want to perform the test suite?"
)

# Diamond download options (only used when running the test suite)
group = parser.add_mutually_exclusive_group(required=False)
group.add_argument(
    "--diamond",
    nargs="?",
    const=DEFAULT_DIAMOND_VERSION,
    default=None,
    help="Force DIAMOND to be downloaded when running the test suite, and optionally specify the version",
)
group.add_argument(
    "--skip-diamond",
    action="store_true",
    default=False,
    help="Skip downloading a DIAMOND binary if running the test suite",
)

args = parser.parse_args()


# ----------------------------------------------------------------------------------------------------------------------#
if __name__ == "__main__":
    # Process input settings
    input_file = args.input_file
    out_file = args.out
    file_type = int(args.file_type)
    org = args.org
    if args.media in ["complete", "rich", "minimal"]:
        media = str(args.media)
    else:
        media = str(args.media).split(",")
    min_frac = float(args.min_frac)
    max_frac = float(args.max_frac)

    if args.tasks is not None:
        metabolic_tasks = str(args.tasks).split(",")
    else:
        metabolic_tasks = None

    new_id = args.name
    if new_id is None:
        new_id = "new_model"

    gram = args.gram
    processors = int(args.cpu) if args.cpu is not None else None
    gapfill = str(args.gapfill)
    exchange_arg = int(args.exchange)
    test = str(args.test)

    # ----------------------------------------------------------------------------------------------------------------------#
    if test == "yes":
        diamond_url = "https://github.com/bbuchfink/diamond/releases"

        # Download a diamond binary if needed
        if args.diamond is not None:
            print(f"Getting DIAMOND v{args.diamond} from {diamond_url}")
            download_diamond(diamond_version=args.diamond)
        if not args.skip_diamond:
            if get_diamond_path() is None:
                print(f"DIAMOND not found...getting DIAMOND from {diamond_url}")
                download_diamond()
            diamond = Diamond()
            print(f"Using DIAMOND v{diamond.get_version()} at {diamond.path}")

        # Download the diamond database file if it hasn't been downloaded yet
        diamond_db_path = resources.get_diamond_db_path()
        if not diamond_db_path.exists():
            resources.download_diamond_db()

        # Run the three tests (each with a different input file)
        # - 488.146.clean.fa: an amino acid .fasta file used to test a type 1 input to reconstructor
        # - JCP8151B.KEGGprot.out: a blast output file used to test a type 2 input to reconstructor
        # - fmt.metaG.01044A.bin.149.KEGGprot.sbml: a .sbml genre used to test a type three input to reconstructor
        test_file_names = [
            "488.146.clean.fa",
            "JCP8151B.KEGGprot.out",
            "fmt.metaG.01044A.bin.149.KEGGprot.sbml",
        ]
        input_types = [1, 2, 3]
        test_num = 0
        for test_file_name, input_type in zip(test_file_names, input_types):
            test_num += 1
            print(f"Performing test {test_num}")

            # Temporary directory to hold test files (and clean them up after test finishes)
            with TemporaryDirectory(dir=resources.RESOURCE_DIR) as tempdir:

                # Extract the test file
                zip_path = resources.RESOURCE_DIR / "testfiles.zip"
                with zipfile.ZipFile(zip_path) as archive:
                    test_file = archive.extract(test_file_name, tempdir)

                # Run reconstructor
                model = reconstruct(
                    test_file, input_type, gram="negative", cpu=processors
                )

                # Save model
                output_file = Path(test_file).with_suffix(".out.sbml")
                cobra.io.write_sbml_model(model, output_file)

        print("Tests finished!")

        quit()
    # ----------------------------------------------------------------------------------------------------------------------#

    # Perform the reconstruction
    model = reconstruct(
        input_file,
        file_type,
        media,
        metabolic_tasks,
        org,
        min_frac,
        max_frac,
        gram,
        new_id,
        processors,
        gapfill == "yes",
        exchange_arg == 1,
    )

    # Write new model to sbml
    if out_file is None:
        input_path = Path(input_file)
        if file_type == 3:
            out_file = input_path.with_suffix(".extended.sbml")
        else:
            out_file = input_path.with_suffix(".sbml")

    print(f"\nSaving new model to {out_file}")
    cobra.io.write_sbml_model(model, out_file)
    print("Finished")
