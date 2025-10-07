"""
This script is used to download the DIAMOND binaries from the DIAMOND releases
page and put them into a single .zip archive in the Reconstructor resources
folder.

This way, the DIAMOND binaries are packaged within reconstructor for all
platforms. Additionally, the .zip archive makes it so they are small enough to
be included in the git repo. Finally, by simply including the .zip archive in
Reconstructor and temporarily extracting the appropriate binary when a user runs
Reconstructor, no files will be left over when Reconstructor gets uninstalled. I
tested the temporary extraction process and it only takes about 0.5s to extract
one of the binaries from the .zip archive, so this doesn't introduce a major
cost in terms of speed.
"""

import tarfile
import zipfile
from tempfile import TemporaryDirectory
import os

import wget

from reconstructor import resources


diamond_release = "v2.1.14"        
url = f"https://github.com/bbuchfink/diamond/releases/download/{diamond_release}/diamond-{{}}"


with TemporaryDirectory() as tempdir:

    # Download and extract Windows exe
    windows = wget.download(url.format("windows.zip"), out=tempdir)
    with zipfile.ZipFile(windows, "r") as zfile:
        zfile.extractall(tempdir)
    windows_exe = os.path.join(tempdir, "diamond-Windows")
    os.rename(os.path.join(tempdir, "diamond.exe"), windows_exe)

    # Download and extract Darwin
    darwin = wget.download(url.format("macos.tar.gz"), out=tempdir)
    with tarfile.open(darwin, "r:gz") as tar:
        tar.extractall(tempdir)
    darwin_exe = os.path.join(tempdir, "diamond-Darwin")
    os.rename(os.path.join(tempdir, "diamond"), darwin_exe)

    #Download and extract Linux
    linux = wget.download(url.format("linux64.tar.gz"), out=tempdir)
    with tarfile.open(linux, "r:gz") as tar:
        tar.extractall(tempdir)
    linux_exe = os.path.join(tempdir, "diamond-Linux")
    os.rename(os.path.join(tempdir, "diamond"), linux_exe)


    # Now place each of the binaries in separate zip archives in the resources folder
    for exe in [windows_exe, darwin_exe, linux_exe]:
        exe_name = os.path.basename(exe)
        zippath = resources.RESOURCE_DIR.joinpath(f"{exe_name}.zip")
        with zipfile.ZipFile(zippath, "w", zipfile.ZIP_DEFLATED) as archive:
            archive.write(exe, exe_name)
