# Reconstructor
# Installation Guide (Must have Python >= 3.8 installed):
## MacOSX
In terminal...
Install Reconstructor package using pip:
```
pip install reconstructor
```
Use following command to run Reconstructor's setup script and final installation checks (should take ~30min to complete depending on your system):
```
python -m reconstructor --test yes
```
## Windows
Install Reconstructor package using CMD.exe prompt launched from Anaconda navigator:
```
pip install reconstructor
```
Use following command to run Reconstructor's setup script and final installation checks (should take ~30min to run depending on your system):
```
python -m reconstructor --test yes
```

# Usage:
## Use reconstructor via COMMAND LINE

* You must run the setup script described in the Installation Guide before proceeding

Now that Reconstructor and all dependency databases are installed, you can proceed to use Reconstructor via command line. An example would be:
```
python -m reconstructor --input_file <input fasta file> --file_type <1,2,3> --gram <negative, positive> --media rich
```
#### Type 1: Build GENRE from annotated amino acid .fasta files
```
python -m reconstructor --input_file Osplanchnicus.aa.fasta --file_type 1 --gram negative --media rich
```

#### Type 2: Build GENRE from BLASTp hits
```
python -m reconstructor --input_file Osplanchnicus.hits.out --file_type 2 --gram negative --media rich
```

#### Type 3: Additional gap-filling (if necessary)
```
python -m reconstructor --input Osplanchnicus.sbml --type 3 --media rich
```
## Use Reconstructor directly in PYTHON
You can use Reconstructor directly in Python for using directly with COBRApy analysis tools.
To import the reconstruction function use the following line: 
``` 
from reconstructor import reconstruct
```
The reconstruct function is defined as follows, with the required argument being input_file. All other arguments have default options following the equals sign. Argument descriptions are provided below
```
reconstruct(input_file, file_type = 1, media=[], org = 'default', min_frac = 0.01, max_frac = 0.5, gram='none', out = 'default', name = 'default', cpu = 1, gapfill = 'yes')
```
Here is an example of how to generate a GENRE from an annotated amino acid fasta file (.fa, type 1 input) directly in your python script.
```
model = reconstruct('218496.4.fa', file_type = 1, gram = 'negative')
```
### Required and optional arguments
```
--input_file <REQUIRED input file, Required, str>
```
```
--file_type <REQUIRED input file type, .fasta = 1, diamond blastp output = 2, .sbml = 3, Required, Default = 1, int> 
```
```
--gram <REQUIRED Type of Gram classificiation (positive or negative), default = positive, str>
```
```
--media <REQUIRED 'rich' is the default and can be used to generate model based on a rich media. List of strings of metabolites in modelseed namespace composing the media condition, comma separated. Must end with _e. For example: 'cpd00001_e'.>
```
```
--org <KEGG organism code. Not required, str>
```
```
--min_frac <Minimum objective fraction required during gapfilling, default = 0.01, float>
```
```
--max_frac <Maximum objective fraction allowed during gapfilling, default = 0.5 ,float>
```
```
--out <Name of output GENRE file, default = default, str>
```
```
--name <ID of output GENRE, default = default, str>
```
```
--cpu <Number of processors to use, default = 1, int>
```

```
--test <run installation tests, default = no, str>
```
