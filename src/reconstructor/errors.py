class ReconstructorError(Exception):
    """
    Base class for Reconstructor errors.
    """

    def __init__(self, *args):
        super().__init__(*args)


class DiamondError(ReconstructorError):
    """
    Base class for errors originating from DIAMOND.
    """
    
    def __init__(self, *args):
        msg = "DIAMOND exited with a non-zero exit status."
        super().__init__(msg, *args)
