# Reconstructor
This repository contains all source code in the reconstructor python package, important file dependencies, and benchmarking scores for reconstructor models. Reconstructor is a COBRApy compatible, automated GENRE building tool from gene fastas based on KEGG annotations.

****Reconstructor is currently only compatible on MacOSX and Linux machines****
#### /RepresentativeGENRES
/RepresentativeGENRES/Reconstructor: Contains 10 representative bacterial GENRES in .sbml format created by reconstructor from annotated .FASTA files.

/RepresentativeGENREs/ModelSEED: Corresponding GENREs created with ModelSEED using the same genome sequence information as used in reconstructor. 

/RepresentativeGRENREs/CarveME: Corresponding GENREs created with CarveME using genome sequences from the same species (or in one case, genus), but not exactly the same strain. These models were taken from the CarveME database ([https://github.com/cdanielmachado/carveme](https://github.com/cdanielmachado/carveme)).

#### /MEMOTE
Contains raw .html files for benchmarking scores for 10 representative reconstructor models

See below in the __Supplementary Model Checks and Analyses__ at the bottom of the README for links to stable .html renderings of the MEMOTE scores and comparisons to other automatic model generation tools

#### /reconstructor
contains all package source code

## Installation:
### Download DIAMOND Aligner
You must first have the diamond sequence aligner downloaded (__MUST BE VERSION 2.0.15__), installation instructions can be found here: [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)  
If you think you already have DIAMOND installed, you can check your version using the command:
```
diamond --version
```

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

# Supplementary Model Checks and Analyses

## Comparison of Reconstructor, ModelSEED, and CarveME GENREs
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

## *C. difficile* reconstructions on Complete (Default), Rich, Minimal, *C. difficile* Defined Minimal, *C. difficile* Defined Enriched

#### /C_difficile_Media
This folder contains .sbml reconstrucitons for *C. difficile* on five different media conditions described below. These include a complete media (default reconstructor media if no media is defined), rich media (built into reconstructor), minimal media (built into reconstructor), *C. difficile* specific defined minimal media (user inputted), and a *C. difficile* specific defined enriched media (user inputted). 

__Complete media__ contains the following ModelSEED compounds:  
cpd00035 (L-Alanine), cpd00051 (L-Arginine), cpd00132 (L-Asparagine), cpd00041 (L-Aspartate), cpd00084 (L-Cysteine), cpd00053 (L-Glutamine), cpd00023 (L-Glutamate), cpd00033 (Glycine), cpd00119 (L-Histidine), cpd00322 (L-Isoleucine), cpd00107 (L-Leucine), cpd00039 (L-Lysine), cpd00060 (L-Methionine), cpd00066 (L-Phenylalanine), cpd00129 (L-Proline), cpd00054 (L-Serine), cpd00161 (L-Threonine), cpd00065 (L-Tryptophan), cpd00069 (L-Tyrosine), cpd00156 (L-Valine), cpd00027 (D-Glucose), cpd00149 (Cobalt), cpd00030 (Manganese), cpd00254 (Magnesium), cpd00971 (Sodium), cpd00063 (Calcium), cpd10515 (Iron), cpd00205 (Potassium), cpd00099 (Choride)  
This reconstruction was generated using the following arugments:  
```
--input 699034.5.fa --type 1 --gram positive  
```
The MEMOTE scores for *C. difficile* reconstruction in complete media can be found [here](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.redo.html)  

__Rich media__ contains the following ModelSEED compounds:  
cpd00001 (Water), cpd00035 (L-Alanine), cpd00041(L-Aspartate), cpd00023 (L-Glutamate), cpd00119 (L-Histidine), cpd00107 (L-Leucine), cpd00060 (L-Methionine), cpd00161 (L-Threonine), cpd00069 (L-Tyrosine), cpd00084 (L-Cysteine), cpd00033 (Glycine), cpd00322 (L-Isoleucine acid), cpd00066 (L-Phenylalanine), cpd00054 (L-Serine), cpd00065 (L-Tryptophan), cpd00156 (L-Valine), cpd00220 (Riboflavin), cpd00644 (Pantothenate), cpd00393 (Folate), cpd00133 (Nicotinamide), cpd00263 (Pyridoxal), cpd00104 (Biotin), cpd00149 (Cobalt), cpd00971 (Sodium), cpd00099 (Chloride), cpd00205 (Potassium), cpd00009 (Phosphate), cpd00063 (Calcium), cpd00254 (Magnesium), cpd10515 (Iron), cpd00030 (Manganese), cpd00242 (Bicarbonate), cpd00226 (Hypoxanthine), cpd01242 (Thyminose), cpd00307 (Cytosine), cpd00092 (Uracil), cpd00117 (D-Alanine), cpd00067 (Hydrogen), cpd00567 (D-Proline), cpd00132 (L-Asparagine), cpd00210 (Taurine), cpd00320 (D-Aspartate), cpd03279 (Deoxyinosine) , cpd00246 (Inosine), cpd00311 (Guanosine), cpd00367 (Cytidine), cpd00277 (Deoxyguanosine), cpd00182 (Adenosine), cpd00654 (Deoxycytidine), cpd00412 (Deoxyuridine), cpd00438 (Deoxyadenosine), cpd00274 (Citrulline), cpd00186 (D-Glutamate), cpd00637 (D-Methionine), cpd00105 (D-Ribose), cpd00305 (Thiamin), cpd00309 (Xanthine), cpd00098 (Choline), cpd00207 (Guanine), cpd00082 (D-Fructose), cpd00129 (L-Proline)
This reconstruction was generated using the following arugments:  
```
--input 699034.5.fa --type 1 --gram positive --media rich_media  
```
The MEMOTE scores for *C. difficile* reconstruction in rich media can be found [here](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.richmedia.html)  

__Minimal media__ contains the following ModelSEED compounds:  
cpd00001 (Water), cpd00065 (L-Tryptophan), cpd00060 (Methionine), cpd00322 (L-Isoleucine), cpd00129 (L-Proline), cpd00156 (L-Valine), cpd00107 (L-Leucine), cpd00084 (L-Cysteine), cpd00149 (Cobalt), cpd00099 (Chloride), cpd10515 (Iron), cpd00030 (Manganese), cpd00254 (Magnesium), cpd00063 (Calcium), cpd00205 (Potassium), cpd00009 (Phosphate), cpd00971 (Sodium), cpd00242 (Carbonate), cpd00104 (Biotin), cpd00644 (Pantothenate), cpd00263 (Pyridoxine) , cpd00082 (D-Fructose) 
This reconstruction was generated using the following arugments:  
```
--input 699034.5.fa --type 1 --gram positive --media minimal_media  
```
The MEMOTE scores for *C. difficile* reconstruction in minimal media can be found [here](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.minimalmedia.html) 

__*C. difficile* defined minimal media__ contains the following ModelSEED compounds:  
cpd00001 (Water), cpd00104 (Biotin), cpd00644 (Pantothenate), cpd00263 (Pyridoxine), cpd00149 (Cobalt), cpd00099 (Chloride), cpd10515 (Iron), cpd00030 (Manganese), cpd00254 (Magnesium), cpd00063 (Calcium), cpd00205 (Potassium), cpd00009 (Phosphate), cpd00971 (Sodium), cpd00242 (Carbonate), cpd00322 (L-Isoleucine), cpd00129 (L-Proline),  cpd00156 (L-Valine), cpd00107 (L-Leucine), cpd00084 (L-Cysteine), cpd00065 (L-Tryptophan), cpd00027 (Glucose)  
This reconstruction was generated using the following arguments:  
```
--input 699034.5.fa --type 1 --gram positive --media 'EX_cpd00001_e','EX_cpd00104_e','EX_cpd00644_e','EX_cpd00263_e','EX_cpd00149_e','EX_cpd00099_e','EX_cpd10515_e','EX_cpd00030_e','EX_cpd00254_e','EX_cpd00063_e','EX_cpd00205_e','EX_cpd00009_e','EX_cpd00971_e','EX_cpd00242_e','EX_cpd00322_e','EX_cpd00129_e','EX_cpd00156_e','EX_cpd00107_e','EX_cpd00084_e','EX_cpd00065_e','EX_cpd00027_e'
```
The MEMOTE scores for *C. difficile* reconstruction in minimal defined media can be found [here](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.definedminimal.html) 

__*C. difficile* defined enriched media__ contains the following ModelSEED compounds:  
cpd00035 (Alanine), cpd00041 (Aspartic Acid), cpd00023 (Glutamic Acid), cpd00119 (Histidine), cpd00107 (Leucine), cpd00060 (Methionine), cpd00129 (Proline), cpd00161 (Threonine), cpd00051 (Arginine), cpd00069 (Tyrosine), cpd00084 (Cysteine), cpd00033 (Glycine), cpd00322 (Isoleucine), cpd00039 (Lysine), cpd00066 (Phenylalanine), cpd00054 (Serine), cpd00065 (Tryptophan), cpd00156 (Valine), cpd00027 (Glucose),
cpd00220 (Riboflavin), cpd00644 (Calcium Pantothentate), cpd00393 (Folate), cpd00133 (Niacin), cpd00263 (Pyridoxine), cpd00104 (Biotin), cpd00149 (Cobalt), cpd00971 (Sodium), cpd00099 (Chloride), cpd00205 (Potassium), cpd00009 (Phosphate), cpd00063 (Calcium),
cpd00254 (Magnesium), cpd10515 (Fe2+), cpd00030 (Mn2+), cpd00242 (Carbonate), cpd00001 (Water), cpd00226 (Hypoxanthine), cpd01242 (Thyminose), cpd00307 (Cytosine), cpd00092 (Uracil)  
This reconstruction was generated using the following arguments" 
```
--input 699034.5.fa --type 1 --gram positive --media 'EX_cpd00035_e','EX_cpd00041_e','EX_cpd00023_e','EX_cpd00119_e','EX_cpd00107_e','EX_cpd00060_e','EX_cpd00129_e','EX_cpd00161_e','EX_cpd00051_e','EX_cpd00069_e','EX_cpd00084_e','EX_cpd00033_e','EX_cpd00322_e','EX_cpd00039_e','EX_cpd00066_e','EX_cpd00054_e','EX_cpd00065_e','EX_cpd00156_e','EX_cpd00027_e','EX_cpd00220_e','EX_cpd00644_e','EX_cpd00393_e','EX_cpd00133_e','EX_cpd00263_e','EX_cpd00104_e','EX_cpd00149_e','EX_cpd00971_e','EX_cpd00099_e','EX_cpd00205_e','EX_cpd00009_e','EX_cpd00063_e','EX_cpd00254_e','EX_cpd10515_e','EX_cpd00030_e','EX_cpd00242_e','EX_cpd00001_e','EX_cpd00226_e','EX_cpd01242_e','EX_cpd00307_e','EX_cpd00092_e'
```
The MEMOTE scores for *C. difficile* reconstruction in enriched defined media can be found [here](https://emmamglass.github.io/ReconstructorMEMOTE.io/699034.5.definedrich.html) 

## Additional Information
Thank you for your interest in reconstructor. If you have any additional questions please email tfz5vy@virginia.edu.

If you encounter any problems, please file an [issue](https://github.com/emmamglass/reconstructor/issues) along with a detailed description.

Distributed under the terms of the [MIT license](https://github.com/emmamglass/reconstructor/blob/main/reconstructor/LICENSE), "reconstructor" is free and open source software
