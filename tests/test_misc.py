"""
This test file is for miscellaneous tests that don't necessarily
belong in any of the other test files.
"""

from copy import deepcopy

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
