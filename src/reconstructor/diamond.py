from typing import Optional, Union, Sequence
import tarfile
import zipfile
from tempfile import TemporaryDirectory
import os
import platform
from urllib.error import HTTPError
import shutil
from importlib import resources
import subprocess
import re

from reconstructor import errors
from reconstructor.utils import download, CallbackT, DownloadProgress


DEFAULT_DIAMOND_VERSION = "2.1.14"
_DIAMOND_URL_TEMPLATE = "https://github.com/bbuchfink/diamond/releases/download/v{version}/diamond-{system}{ext}"
_WINDOWS_BIN = "diamond.exe"
_MACOS_BIN = "diamond"
_LINUX_BIN = "diamond"

_BIN_DIR = resources.files(__package__).joinpath("bin")


class Diamond:

    def __init__(self, bin_path: Optional[Union[str, bytes, os.PathLike]] = None):
        if bin_path is None:
            bin_path = get_diamond_path()
            if bin_path is None:
                raise errors.DiamondNotFoundError()
        self.path: str = str(bin_path)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.path)})"

    def __str__(self) -> str:
        return f"DIAMOND v{self.get_version()} ({self.path})"

    def __call__(
        self, options: Sequence[str], *args, **kwargs
    ) -> subprocess.CompletedProcess:
        options = [self.path] + options
        result: subprocess.CompletedProcess = subprocess.run(options, *args, **kwargs)
        try:
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            raise errors.DiamondProcessError(e) from e
        else:
            return result

    def blastp(
        self,
        db: Union[str, bytes, os.PathLike],
        query: Union[str, bytes, os.PathLike],
        out: Optional[Union[str, bytes, os.PathLike]] = None,
        *options: str,
        **sp_kwargs,
    ):
        args = ["blastp", "--db", db, "--query", query]
        if out is not None:
            args.extend(["--out", out])
        args.extend([str(x) for x in options])
        return self.__call__(args, **sp_kwargs)

    def get_version(self):
        args = ["version"]
        result: subprocess.CompletedProcess[str] = self.__call__(
            args, capture_output=True, text=True
        )
        pattern = re.compile(r"diamond version (\d+\.\d+\.\d+)")
        version = pattern.sub(r"\1", result.stdout.strip())
        return version


def get_diamond_path(name: Optional[str] = None) -> Optional[str]:
    """
    Get the path to a DIAMOND executable or return None if one is not found.

    Checks the following locations (in the order shown) and returns the first
    one that is found:

    1. reconstructor/bin/ (a bin folder inside the reconstructor package directory)
    2. The user's PATH by calling `shutil.which(name)`

    If DIAMOND is not found, then `None` is returned.
    """
    if name is None:
        name = "diamond.exe" if platform.system() == "Windows" else "diamond"

    # Check the bin folder inside reconstructor
    path = _BIN_DIR.joinpath(name)
    if os.path.exists(path):
        return str(path)

    # Get from PATH
    return shutil.which(name)


def download_diamond(
    dir: Optional[Union[str, bytes, os.PathLike]] = _BIN_DIR,
    diamond_version: str = DEFAULT_DIAMOND_VERSION,
    bin_name: Optional[str] = None,
    callback: Optional[CallbackT] = DownloadProgress("Downloading DIAMOND..."),
) -> str:
    """
    Download the appropriate DIAMOND binary for the operating system from
    [GitHub](https://github.com/bbuchfink/diamond/releases).
    """
    if dir is _BIN_DIR and not os.path.exists(_BIN_DIR):
        os.mkdir(_BIN_DIR)

    system = platform.system()
    if system == "Windows":
        return _download_windows(dir, diamond_version, bin_name, callback)

    if dir is None:
        dir = ""

    if system == "Darwin":
        return _download_macos(dir, diamond_version, bin_name, callback)
    if system == "Linux":
        return _download_linux(dir, diamond_version, bin_name, callback)
    # TODO: is it possible for system to be something else?


def cleanup_bin():
    """
    Deletes the bin directory inside the reconstructor package.
    """
    if os.path.exists(_BIN_DIR):
        shutil.rmtree(_BIN_DIR)


def _download_windows(
    dir: Optional[Union[str, bytes, os.PathLike]],
    diamond_version: str,
    bin_name: Optional[str],
    callback: Optional[CallbackT] = None,
) -> str:
    if bin_name is None:
        bin_name = _WINDOWS_BIN

    with TemporaryDirectory() as tempdir:
        download_path = _download_archive(
            tempdir, diamond_version, "windows", ".zip", callback
        )
        with zipfile.ZipFile(download_path, "r") as zfile:
            zfile.extract(bin_name, dir)

    return os.path.join(dir, bin_name)


def _download_macos(
    dir: Optional[Union[str, bytes, os.PathLike]],
    diamond_version: str,
    bin_name: Optional[str],
    callback: Optional[CallbackT] = None,
) -> str:
    if bin_name is None:
        bin_name = _MACOS_BIN

    with TemporaryDirectory() as tempdir:
        download_path = _download_archive(
            tempdir, diamond_version, "macos", ".tar.gz", callback
        )
        with tarfile.open(download_path, "r:gz") as tar:
            tar.extract(bin_name, dir)

    return os.path.join(dir, bin_name)


def _download_linux(
    dir: Optional[Union[str, bytes, os.PathLike]],
    diamond_version: str,
    bin_name: Optional[str],
    callback: Optional[CallbackT] = None,
) -> str:
    if bin_name is None:
        bin_name = _LINUX_BIN

    with TemporaryDirectory() as tempdir:
        download_path = _download_archive(
            tempdir, diamond_version, "linux64", ".tar.gz", callback
        )
        with tarfile.open(download_path, "r:gz") as tar:
            tar.extract(bin_name, dir)

    return os.path.join(dir, bin_name)


def _download_archive(
    dir: Optional[Union[str, bytes, os.PathLike]],
    diamond_version: str,
    system: str,
    ext: str,
    callback: Optional[CallbackT] = None,
) -> str:
    url = _DIAMOND_URL_TEMPLATE.format(version=diamond_version, system=system, ext=ext)
    download_path = os.path.join(dir, "diamond.zip")

    try:
        download(url, download_path, callback)
    except HTTPError as e:
        raise errors.DiamondDownloadError(e, diamond_version, system) from e

    return download_path
