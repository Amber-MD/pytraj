import os
from contextlib import contextmanager
import tempfile
from shutil import rmtree

try:
    from ..externals.wurlitzer import pipes
except ImportError:

    def pipes():
        # win sucks
        yield "", ""


capture_stdout = pipes


@contextmanager
def tempfolder():
    my_temp = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(my_temp)
    yield
    os.chdir(cwd)
    rmtree(my_temp)
