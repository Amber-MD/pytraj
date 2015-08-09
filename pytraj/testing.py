from __future__ import absolute_import, print_function
import os

from .decorators import test_if_having, no_test, local_test
from .datafiles.load_sample_data import load_sample_data
from .utils import eq, aa_eq
from .utils.check_and_assert import is_linux
from .utils import duplicate_traj, Timer
from .Topology import Topology
from .utils.context import goto_temp_folder

__all__ = ['test_if_having', 'no_test', 'local_test', 'load_sample_data', 'eq',
           'aa_eq', 'is_linux', 'make_random_frame', 'duplicate_traj',
           'make_fake_traj', 'Timer', 'amberhome', 'cpptraj_test_dir',
           'run_docstring']

try:
    amberhome = os.environ['AMBERHOME']
    cpptraj_test_dir = os.path.join(amberhome, 'AmberTools', 'test', 'cpptraj')
except:
    amberhome = None
    cpptraj_test_dir = None

possible_path = "../cpptraj/test/"
if os.path.exists(possible_path):
    cpptraj_test_dir = possible_path


def make_random_frame(n_atoms=10000):
    import numpy as np
    from pytraj import Frame

    frame = Frame(n_atoms)
    frame.xyz[:] = np.random.randn(n_atoms, 3)
    return frame


header_doc = '''
from pytraj import io
import pytraj.common_actions as pyca
traj = io.load_sample_data("tz2")
'''


def run_docstring(func):
    func_doct = func.__doc__
    _doc = [x.lstrip() for x in func.__doc__.split("\n")]
    _doc = filter(lambda x: '>>>' in x, _doc)
    _doc = [x.replace(">>> ", "") for x in _doc]
    doc = "\n".join(x for x in _doc)
    if "from pytraj import io" not in func_doct:
        doc = "\n".join((header_doc, doc))
    exec(doc)


if __name__ == "__main__":
    print(amberhome)
    print(cpptraj_test_dir)
