from __future__ import absolute_import, print_function
import os

from .decorators import test_if_having, no_test

try:
    amberhome = os.environ['AMBERHOME']
    cpptraj_test_dir = os.path.join(amberhome, 'AmberTools', 'test', 'cpptraj')
except:
    amberhome = None
    cpptraj_test_dir = None

if __name__ == "__main__":
    print (amberhome)
    print (cpptraj_test_dir)
