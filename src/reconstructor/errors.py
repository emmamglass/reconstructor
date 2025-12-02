import subprocess
from urllib.error import HTTPError


class ReconstructorError(Exception):
    """
    Base class for Reconstructor errors.
    """

    pass


class DiamondError(ReconstructorError):
    """
    Base class for errors originating from DIAMOND.
    """

    pass


class DiamondNotFoundError(DiamondError, FileNotFoundError):
    """
    The error raised if a DIAMOND executable cannot be found.
    """

    def __init__(self):
        msg = (
            "DIAMOND was not found. If you have already installed it, make sure it is available on "
            "your PATH. Otherwise you can install it by running the DIAMOND test suite:\n\n\tpython "
            "-m reconstructor --test yes\n"
        )
        super().__init__(msg)


class DiamondDownloadError(DiamondError, HTTPError):
    """
    The error raised if an HTTPError is encountered while trying to download
    DIAMOND.
    """

    def __init__(self, error: HTTPError, version, system):
        super().__init__(error.url, error.code, error.msg, error.hdrs, error.fp)
        self.version = version
        self.system = system

    def __str__(self) -> str:
        return (
            f"DIAMOND v{self.version} for {self.system} was not found. Please double check that a binary "
            "for this version and system is available at https://github.com/bbuchfink/diamond/releases."
        )


class DiamondProcessError(DiamondError, subprocess.CalledProcessError):
    """
    The error that gets raised if DIAMOND encounters a CalledProcessError.
    """

    def __init__(self, error: subprocess.CalledProcessError):
        super().__init__(
            returncode=error.returncode,
            cmd=error.cmd,
            output=error.output,
            stderr=error.stderr,
        )
