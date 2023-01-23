# Reconstructor
This repository contains all source code in the reconstructor python package, important file dependencies, and benchmarking scores for reconstructor models. Reconstructor is a COBRApy compatible, automated GENRE building tool from gene fastas based on KEGG annotations.

****Reconstructor is currently only compatible on MacOSX and Linux machines****
#### /RepresentativeGENRES
/RepresentativeGENRES/Reconstructor: Contains 10 representative bacterial GENRES in .sbml format created by reconstructor from annotated .FASTA files.

/RepresentativeGENREs/ModelSEED: Corresponding GENREs created with ModelSEED using the same genome sequence information as used in reconstructor. 

/RepresentativeGRENREs/CarveME: Corresponding GENREs created with CarveME using genome sequences from the same species (or in one case, genus), but not exactly the same strain. These models were taken from the CarveME database ([https://github.com/cdanielmachado/carveme](https://github.com/cdanielmachado/carveme)).

#### /MEMOTE
Contains raw .html files for benchmarking scores for 10 representative reconstructor models

To view rendered .html MEMOTE benchmarking scores for Reconstructor models and the corresponding ModelSEED/CarveME models, use the links below:  
__*Porphyromonas gingivalis:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/837.134.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/837.83MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/837.83C.html)  
__*Bacillus amyloliquefaciens:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/1390.556.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/1390.556MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/1390.556C.html)  
__*Citrobacter braakii:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/57706.84.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/57706.84MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/57706.84C.html)  
__*Acinetobacter bereziniae:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/106648.24.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/106648.24MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/106648.24C.html)  
__*Tropheryma whipplei:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/218496.4.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/218496.4MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/218496.4C.html)  
__*Proteus mirabilis:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/529507.6.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/529507.6MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/529507.6C.html)  
__*Clostridiodes difficile:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5C.html)  
__*Campylobacter jejuni:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/1349827.3.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/1349827.3MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/1349827.3C.html)  
__*Helicobacter pylori:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/1382925.3.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/1382925.3MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/1382925.3C.html)  
__*Escherichia coli:*__ [Resonstructor](https://emmamglass.github.io/ReconstructorMEMOTE.io/2848143.3.redo.html) || [ModelSEED](https://emmamglass.github.io/ReconstructorMEMOTE.io/2848143.3MS.html) || [CarveME](https://emmamglass.github.io/ReconstructorMEMOTE.io/2848143.3C.html)  

Reconstructor outperforms ModelSEED and CarveMe in terms of memote scores in all MEMOTE subcategories and in overall total score:
<p align="center">
  <img src="https://github.com/emmamglass/reconstructor/blob/main/MemoteComparison.png" width=70% height=70% class="center">
</p>

A) MEMOTE consistency score break down by category, B) MEMOTE score for each overall category, C) Total overall score comparison

#### /reconstructor
contains all package source code

## Installation:
### Download DIAMOND Aligner
You must first have the diamond sequence aligner downloaded (__MUST BE VERSION 2.0.15__), installation instructions can be found here: [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)

### Install Reconstructor python package
This can be done via pip in terminal on mac

```
pip install reconstructor
```

*You must be running >= Python 3.8

Some windows configurations may require you to include ``` --user ``` in the pip install line like below:
```
pip install --user reconstructor
```

## Test suite:
#### Use the following command to run the test suite
Run the following test to ensure reconstruction was installed correctly and is functional. This series of tests should take about an hour to run, dependent on computer/processor speed. These are runtimes for reconstructor on a 2020 MacBook Pro with a 1.4 GHz Quad-Core Intel Core i5 processor.

MAC USERS MAY BE ASKED FOR TERMINAL TO HAVE ACCESS TO DOWNLOADS, CAMERA, LOCATION, ETC. Please allow terminal to have access to all locations on your computer. Reconstructor will NOT gather data from your camera, location, or other sensitive information. Reconstructor is simply searching for the file titled glpk_interface.py on your local machine (installed when COBRA module is installed) and replacing it with a newer, functional version.

Use the command below to test reconstructor to ensure correct installation. 

```
python -m reconstructor --test yes
```
*YOU MUST RUN THE TEST SUITE BEFORE PROCEEDING TO USE RECONSTRUCTOR

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
### Required and optional arguments
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

```
--test <run installation tests, default = no>
```
## Additional Information
Thank you for your interest in reconstructor. If you have any additional questions please email tfz5vy@virginia.edu.

If you encounter any problems, please file an [issue](https://github.com/emmamglass/reconstructor/issues) along with a detailed description.

Distributed under the terms of the [MIT license](https://github.com/emmamglass/reconstructor/blob/main/reconstructor/LICENSE), "reconstructor" is free and open source software
