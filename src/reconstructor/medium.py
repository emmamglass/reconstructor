from typing import Sequence


_media_dict: dict[str, list[str]] = {}


def get_medium(name: str):
    global _media_dict
    medium = _media_dict.get(name)
    if medium is None:
        raise KeyError(f"{repr(name)} is not a registered medium.")
    return medium


def register(medium: Sequence[str], name: str):
    global _media_dict
    _media_dict[name] = list(medium)


RICH = [
    "cpd00001_e",
    "cpd00035_e",
    "cpd00041_e",
    "cpd00023_e",
    "cpd00119_e",
    "cpd00107_e",
    "cpd00060_e",
    "cpd00161_e",
    "cpd00069_e",
    "cpd00084_e",
    "cpd00033_e",
    "cpd00322_e",
    "cpd00066_e",
    "cpd00054_e",
    "cpd00065_e",
    "cpd00156_e",
    "cpd00220_e",
    "cpd00644_e",
    "cpd00393_e",
    "cpd00133_e",
    "cpd00263_e",
    "cpd00104_e",
    "cpd00149_e",
    "cpd00971_e",
    "cpd00099_e",
    "cpd00205_e",
    "cpd00009_e",
    "cpd00063_e",
    "cpd00254_e",
    "cpd10515_e",
    "cpd00030_e",
    "cpd00242_e",
    "cpd00226_e",
    "cpd01242_e",
    "cpd00307_e",
    "cpd00092_e",
    "cpd00117_e",
    "cpd00067_e",
    "cpd00567_e",
    "cpd00132_e",
    "cpd00210_e",
    "cpd00320_e",
    "cpd03279_e",
    "cpd00246_e",
    "cpd00311_e",
    "cpd00367_e",
    "cpd00277_e",
    "cpd00182_e",
    "cpd00654_e",
    "cpd00412_e",
    "cpd00438_e",
    "cpd00274_e",
    "cpd00186_e",
    "cpd00637_e",
    "cpd00105_e",
    "cpd00305_e",
    "cpd00309_e",
    "cpd00098_e",
    "cpd00207_e",
    "cpd00082_e",
    "cpd00129_e",
]

COMPLETE = [
    "cpd00035_e",
    "cpd00051_e",
    "cpd00132_e",
    "cpd00041_e",
    "cpd00084_e",
    "cpd00053_e",
    "cpd00023_e",
    "cpd00033_e",
    "cpd00119_e",
    "cpd00322_e",
    "cpd00107_e",
    "cpd00039_e",
    "cpd00060_e",
    "cpd00066_e",
    "cpd00129_e",
    "cpd00054_e",
    "cpd00161_e",
    "cpd00065_e",
    "cpd00069_e",
    "cpd00156_e",
    "cpd00027_e",
    "cpd00149_e",
    "cpd00030_e",
    "cpd00254_e",
    "cpd00971_e",
    "cpd00063_e",
    "cpd10515_e",
    "cpd00205_e",
    "cpd00099_e",
]

MINIMAL = [
    "cpd00001_e",
    "cpd00065_e",
    "cpd00060_e",
    "cpd00322_e",
    "cpd00129_e",
    "cpd00156_e",
    "cpd00107_e",
    "cpd00084_e",
    "cpd00149_e",
    "cpd00099_e",
    "cpd10515_e",
    "cpd00030_e",
    "cpd00254_e",
    "cpd00063_e",
    "cpd00205_e",
    "cpd00009_e",
    "cpd00971_e",
    "cpd00242_e",
    "cpd00104_e",
    "cpd00644_e",
    "cpd00263_e",
    "cpd00082_e",
]


# Register the media
register(RICH, "rich")
register(COMPLETE, "complete")
register(MINIMAL, "minimal")
