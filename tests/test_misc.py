"""
This test file is for miscellaneous tests that don't necessarily
belong in any of the other test files.
"""

from copy import deepcopy
import itertools

import pytest
import cobra

from . import utils


def test_time_to_copy_reaction(universal_model: cobra.Model):
    """
    This tests how long it takes to make a deepcopy of a reaction in the
    universal model.
    
    Since the model is currently stored as a pickle file, it is more sensitive
    to changes in the version of cobrapy. For example, using cobra v0.29 makes
    it take 30+ seconds to copy a single reaction. This can make running
    Reconstructor prohibitevly slow, so this test is to ensure that 
    """
    timed_copy = utils.add_timer(deepcopy, timeout=5)   # Timeout to keep the test
                                                        # test from running for a
                                                        # really long time
    
    n_replicates = 5
    times = []
    for i in range(n_replicates):
        rxn, time = timed_copy(universal_model.reactions[0]) 
        times.append(time)
    
    avg_time = sum(times) / len(times)
    assert avg_time < 0.1


def test_write_sbml(universal_model: cobra.Model, tmp_path):
    """
    Tests writing a model to SBML (the final action of the Reconstructor script)
    that has reactions copied from the universal model.
    
    With newer versions of cobrapy, reactions from the universal model are
    sometimes missing attributes that cause errors when writing the model to
    SBML. The `_annotation` is added in the Reconstructor script so we add it
    here as well.
    """
    new_model = cobra.Model("test_model")

    # Copy a reaction from the universal model
    new_model.add_reactions([deepcopy(universal_model.reactions[0])])

    # Add _annotation to all genes, reactions, and metabolites (causes error otherwise)
    for obj in itertools.chain(new_model.genes, new_model.metabolites, new_model.reactions):
        obj._annotation = {}

    # Test writing the model
    model_path = tmp_path / "model.sbml"
    cobra.io.write_sbml_model(new_model, str(model_path))
