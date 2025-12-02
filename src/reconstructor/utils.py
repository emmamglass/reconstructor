from typing import Callable, Union, Optional, Any
import os
import re
from urllib import request
import http.client
from tempfile import TemporaryDirectory
import time


CallbackT = Callable[[int, int, int], Any]


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


def download(
    url: str, path: Union[str, bytes, os.PathLike], callback: Optional[CallbackT] = None
):
    """
    Download the contents of a url and save to the specified path.
    """

    with TemporaryDirectory() as tmpdir:
        tmp_path = os.path.join(tmpdir, "download.tmp")

        # Download contents to a temporary path
        with request.urlopen(url) as response, open(tmp_path, "wb") as file:
            response: http.client.HTTPResponse
            total_size = int(response.info().get("content-length", 0))
            block_size = 64 * 1024
            count = 0

            while True:
                chunk = response.read(block_size)
                if not chunk:
                    break

                count += 1
                file.write(chunk)

                if callback is not None:
                    callback(count, block_size, total_size)

        # Rename the temporary downloaded file to the provided path
        os.replace(tmp_path, path)

    return path


class DownloadProgress:
    """
    A callback object to display the progress of a download.
    """

    def __init__(self, msg: str = "Downloading...", freq: float = 0.05):
        self.msg: str = msg
        self.freq: float = freq
        self._prev: Optional[float] = None

    def __call__(self, count: int, block: int, total: int) -> None:
        now = time.monotonic()
        finished = count * block >= total
        if self._prev is None or finished or (now - self._prev) >= self.freq:
            self._prev = now
            progress = f"\r{self.msg} {count*block/total:.1%}"
            end = "\n" if finished else ""
            print(progress, end=end, flush=True)
