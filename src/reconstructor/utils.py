import re


def is_valid_sbml_id(id_: str) -> bool:
    """
    Checks if an ID is a valid SBML ID.

    SBML IDs should start either with an underscore (`_`) or a letter and should
    only contain underscores, letters, and numbers.
    """
    sbml_pattern = r"^[_a-zA-Z]\w*$"
    return bool(re.match(sbml_pattern, id_))


def sanitize_sbml_id(id_: str):
    """
    Formats an ID to be a valid SBML ID.

    Replaces any invalid characters (e.g., spaces, dashes, etc.) with `_` and
    adds `_` at the beginning if the string starts with a number. Note that if a
    string contains a sequence of multiple invalid characters in a row, the
    sequence of characters will be replaced by a single underscore.

    Examples
    --------
    >>> sanitize_sbml_id("a_valid_id")
    'a_valid_id'
    >>> sanitize_sbml_id("an invalid--id #3")
    'an_invalid_id_3'
    >>> sanitize_sbml_id("3-atp")
    '_3_atp'
    """
    if len(id_) > 0 and id_[0].isdigit():
        id_ = "_" + id_
    invalid_chars = r"\W+"
    return re.sub(invalid_chars, "_", id_)
