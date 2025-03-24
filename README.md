# Reconstructor - Current Version: 1.1.2

This repository contains all source code in the `Reconstructor` Python package,
important file dependencies, and benchmarking scores for `Reconstructor`
models. `Reconstructor` is a COBRApy compatible, automated GENRE building tool
from annotated aminoa acid .fasta files based on KEGG annotations.

> [!NOTE]
> `Reconstructor` is currently only compatible on MacOSX and Windows machines

## Quick Installation Guide

### MacOSX pre-installation steps

In a terminal....

Install homebrew if you do not already have it installed by copying/pasting the
following line into the terminal (further instructions [here](https://brew.sh)):

```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Download the DIAMOND sequence aligner (version >= 2.0.15) using homebrew:

```shell
brew install diamond
```

Now continue to the [instructions for all platforms](#installation-all-platforms).

### Installation (all platforms)

Install `Reconstructor` using pip:

```shell
pip install reconstructor
```

Use following command to download additional dependencies/run final installation
checks (should take ~1hr to complete depending on your system):

```shell
python -m reconstructor --test yes
```

See the [documentation](docs/installation.md) for additional details about
installation.

## Usage

### Use reconstructor via COMMAND LINE

Now that `Reconstructor` and all dependency databases are installed, you can
proceed to use `Reconstructor` via the command line. An example would be:

```shell
python -m reconstructor --input_file <input fasta file> --file_type <1,2,3> --gram <negative, positive> --media rich
```

All possible command line arguments are described in a
[later section](#required-and-optional-arguments).

#### Type 1: Build GENRE from annotated amino acid .fasta file

```shell
python -m reconstructor --input_file Osplanchnicus.aa.fasta --file_type 1 --gram negative --media rich
```

#### Type 2: Build GENRE from BLASTp hits

```shell
python -m reconstructor --input_file Osplanchnicus.hits.out --file_type 2 --gram negative --media rich
```

#### Type 3: Additional gap-filling (if necessary)

```shell
python -m reconstructor --input Osplanchnicus.sbml --type 3 --media rich
```

### Use Reconstructor directly in PYTHON

Reconstructor can be imported and used directly in Python for easy integration
with COBRApy analysis tools:

```python
from reconstructor import reconstruct
```

Here is an example of how to generate a GENRE from an annotated amino acid fasta
file (.fa, type 1 input) directly in your python script.

```python
model = reconstruct('218496.4.fa', file_type=1, gram='negative')
```

All arguments for this function are described in the
[next section](#required-and-optional-arguments).

### Required and optional arguments

```shell
--input_file <REQUIRED input file, required, str>
```

```shell
--file_type <REQUIRED input file type, .fasta = 1, diamond blastp output = 2, .sbml = 3, Required, Default = 1, int>
```

```shell
--gram <REQUIRED Type of Gram classificiation (positive or negative), default = positive, str>`
```

```shell
--media <REQUIRED 'rich' is the default and can be used to generate model based on a rich media. List of strings of metabolites in modelseed namespace composing the media condition, comma separated. Must end with _e. For example: 'cpd00001_e'.>`
```

```shell
--org <KEGG organism code. Not required, str>
```

```shell
--min_frac <Minimum objective fraction required during gapfilling, default = 0.01, float>
```

```shell
--max_frac <Maximum objective fraction allowed during gapfilling, default = 0.5 ,float>
```

```shell
--out <Name of output GENRE file, default = default, str>
```

```shell
--name <ID of output GENRE, default = default, str>
```

```shell
--cpu <Number of processors to use, default = 1, int>
```

```shell
--test <run installation tests, default = no, str>
```

## Git Repository Structure

### [reconstructor](reconstructor/)

Contains the source code for `Reconstructor`.

### [reference](reference/)

Contains a variety of zipped reference files related to `Receonstructor`:

- RepresentativeGENRES.zip (comparison of models created by `Reconstructor`,
  ModelSeed, and CarveME)
- MEMOTE.zip (MEMOTE benchmarking scores for 10 representative `Reconstructor`
  models)
- C_difficile_media.zip (*C. Difficile* models made with 5 different media)

See the [reference folder README](reference/README.md) for more details about
its contents.

### [docs](docs/)

Contains documentation files for `Reconstructor`.

> [!WARNING]
> The documentation is still very much under development.

## Additional Information

Thank you for your interest in `Reconstructor`. If you have any additional
questions please open an issue.

If you encounter any problems, please file an [issue][github_issues] along with
a detailed description.

Distributed under the terms of the [MIT license](LICENSE), `Reconstructor` is
free and open source software.

[github_issues]: https://github.com/emmamglass/reconstructor/issues
