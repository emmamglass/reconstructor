# Reconstructor
Reconstructor is a COBRApy compatible, automated GENRE building tool from gene fastas based on KEGG annotations. For necessary installation files go to: https://github.com/emmamglass/reconstructor/

#### /MEMOTE (on github)
contains benchmarking scores for 10 representative reconstructor models

#### /reconstructor (on github)
contains all package source code

#### /testfiles (on github)
folder containing files necessary to test for successful installation

## Installation:
### 1) Install Reconstructor python package
This can be done via pip in the command line

```
pip install reconstructor
```

### 2) Download necessary reference databases
Go to https://github.com/emmamglass/reconstructor/releases/tag/v0.0.1 and download all assets (excluding source code zip files)

### 3) Create 'refs' folder
Go to your local downloads folder and create a folder called 'refs' containing the downloaded assets:
```
biomass.sbml 
```
```
compounds.json
```
```
gene_modelseed.pickle
```
```
gene_names.pickle
```
```
reactions.json
```
```
screened_kegg_prokaryotes_pep_db.dmnd
```
```
universal.pickle
```

### 4) scp refs folder into the reconstructor package folder
Use the following command (or similar) in mac terminal to copy the refs folder into the reconstructor python package folder
```
scp -r ~/Downloads/refs ~/opt/anaconda3/lib/python3.9/site-packages/reconstructor
```

### 5) Download diamond 
Diamond version v2.0.15 or higher is REQUIRED. Install instructions for diamond can be found here if compiling from source: https://github.com/bbuchfink/diamond/wiki/2.-Installation. 

Alternatively, diamond can be installed via homebrew:
https://formulae.brew.sh/formula/diamond

Diamond must be v2.0.15 or higher.

## Test suite:
#### 1) Download testfiles folder to your desktop
Download the testfiles folder onto your local desktop.

#### 2) Change directory to testfiles folder
Use command line to change your current directory to the testfiles folder. Use the following command or something similar:
```
cd Desktop/testfiles
```

#### 3) Run the following tests to ensure installation was successful
Run the following three tests to ensure reconstruction was installed correctly and is functional. The first test will take approximately 45 minutes to run, second test ~8  minutes, third test ~2 minutes, dependent on computer/processor speed. :
#### Test 1
```
python -m reconstructor --input 488.146.fa --type 1 --gram negative
```
#### Test 2
```
python -m reconstructor --input 488.146a.KEGGprot.out --type 2 --gram negative
```
#### Test 3
```
python -m reconstructor --input fmt.metaG.01044A.bin.149.KEGGprot.sbml --type 3
```

#### 4) Delete testfiles folder from your destktop
You no longer need the testfiles folder after running the installation tests, so feel free to delete this from your desktop.

## Usage:
### Use reconstructor via command line
Now that reconstructor and all dependency databases are installed, you can proceed to use reconstructor via command line. An example would be:
```
python -m reconstructor --input <input fasta file> --type <1,2,3> --gram <negative, positive> --other arguments <args>
```
#### Type 1: Build GENRE from annotated amino acid fasta files
```
python -m reconstructor --input Osplanchnicus.aa.fasta --type 1 --gram negative --other_args <args>
```

#### Type 2: Build GENRE from BLASTp hits
```
python -m reconstructor --input Osplanchnicus.hits.out --type 2 --gram negative --other_args <args>
```

#### Type 3: Additional gap-filling (if necessary)
```
python -m reconstructor --input Osplanchnicus.sbml --type 3 --other_args <args>
```
### Required and optional parameters
```
--input <input file, Required>
```
```
--type <input file type, .fasta = 1, diamond blastp output = 2, .sbml = 3, Required, Default = 1> 
```
```
--gram <Type of Gram classificiation (positive or negative), default = positive>
```
```
--media <List of metabolites composing the media condition. Not required.>
```
```
--tasks <List of metabolic tasks. Not required>
```
```
--org <KEGG organism code. Not required>
```
```
--min_frac <Minimum objective fraction required during gapfilling, default = 0.01>
```
```
--max_frac <Maximum objective fraction allowed during gapfilling, default = 0.5>
```
```
--out <Name of output GENRE file, default = default>
```
```
--name <ID of output GENRE, default = default>
```
```
--cpu <Number of processors to use, default = 1>
```