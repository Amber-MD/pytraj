import os
from contextlib import contextmanager
import tempfile
from shutil import rmtree


@contextmanager
def goto_temp_folder():
    my_temp = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(my_temp)
    yield
    os.chdir(cwd)
    rmtree(my_temp)
