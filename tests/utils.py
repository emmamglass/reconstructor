import time
import functools
from typing import Optional
import signal
import platform
from contextlib import contextmanager


@contextmanager
def _timeout(duration: float):
     def timeout_handler(signum, frame):
          raise TimeoutError(f"Timed out after {duration} seconds")
     
     signal.signal(signal.SIGALRM, timeout_handler)
     signal.alarm(duration)
     try:
          yield
     finally:
          signal.alarm(0)


def add_timer(f, timeout: Optional[float] = None):
    """
    Creates a wrapper that will return the result of the wrapped function
    along with the elapsed time in seconds.

    Additionally, a timeout can be set to stop the execution of the
    wrapped function after the specified number of seconds. Note that the
    timeout functionality will not work on Windows (becuase it uses
    signal.SIGALRM). Consequently, specified timeouts will be ignored
    if the platform is detected to be Windows.
    """

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
            if platform.system() == "Windows" or timeout is None:
                start = time.time()
                result = f(*args, **kwargs)

            else:
                 with _timeout(timeout):
                      start = time.time()
                      result = f(*args, **kwargs)
                      
            elapsed = time.time() - start
            return result, elapsed
    
    return wrapper
