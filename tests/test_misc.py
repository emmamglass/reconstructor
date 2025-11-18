"""
This test file is for miscellaneous tests that don't necessarily
belong in any of the other test files.
"""


import pytest
import cobra

from reconstructor import medium


def test_write_sbml(universal_model: cobra.Model, tmp_path):
    """
    Tests writing a model to SBML (the final action of the Reconstructor script)
    that has reactions copied from the universal model.
    
    With newer versions of cobrapy, reactions from the universal model are
    sometimes missing attributes that cause errors when writing the model to
    SBML.
    """
    new_model = cobra.Model("test_model")

    # Copy a reaction from the universal model
    new_model.add_reactions([universal_model.reactions[0].copy()])

    # Test writing the model
    model_path = tmp_path / "model.sbml"
    cobra.io.write_sbml_model(new_model, str(model_path))


def test_medium():
    test_medium = ["a", "b", "c", "d", "e"]
    medium.register(test_medium, "test")
    assert test_medium == medium.get_medium("test")


def test_medium_error():
    with pytest.raises(KeyError):
        medium.get_medium("missing")
