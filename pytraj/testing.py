from __future__ import absolute_import, print_function
import os

from .decorators import test_if_having, no_test
from .data_sample.load_sample_data import load_sample_data
from .utils import eq, aa_eq

try:
    amberhome = os.environ['AMBERHOME']
    cpptraj_test_dir = os.path.join(amberhome, 'AmberTools', 'test', 'cpptraj')
except:
    amberhome = None
    cpptraj_test_dir = None

if __name__ == "__main__":
    print (amberhome)
    print (cpptraj_test_dir)
