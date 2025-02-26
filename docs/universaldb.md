# User-Specified Universal Database Modification

The universal reaction database used in `Reconstructor` is a modified version of
the [ModelSEED biochemistry database](https://github.com/ModelSEED/ModelSEEDDatabase).
We generated the universal database using the
[UniversalCreation.py](src/UniversalCreation.py) script. We corrected poorly
defined metabolite formulas and removed mass imbalanced reactions using the
[curateuniversal.py](src/curateuniversal.py) script that is avaliable in this
repository.  

If you wish to modify the provided universal reaction database you must first
locate the universal.pickle file that was downloaded during the `Reconstructor`
installation phase on your local computer. You can modify this file by running a
python script. You must write this script yourself or use
[curateuniversal.py](src/curateuniversal.py) as a template. You can then run
this script in your command line.

## Dependencies and loading universal.pickle

Begin your script with the following dependencies:

```python
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
```

Then, load your universal.pickle file using:

```python
with open("path/to/universal.pickle", "rb") as f:
    universal = pickle.load(f)
```

## Maniuplating the universal database

The loaded `universal` variable is a COBRApy `Model` object which can now be
manipulated. This means you can add, change, or remove reactions and metabolites
from the universal reaction database in the same way that you would with any
other GENRE in COBRApy. For information about working with COBRApy models,
please refer to the [COBRApy documentation][cobrapy_docs].

[cobrapy_docs]: https://cobrapy.readthedocs.io/en/latest/building_model.html  

## Saving the new universal database

After making modifications to the universal database, you can replace the
existing universal database using the following command or similar:

```python
with open("path/to/universal.pickle", "wb") as f:
    pickle.dump(universal, f)
```
