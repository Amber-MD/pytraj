import os

try:
    amberhome = os.environ['AMBERHOME']
    has_amberhome = True
    cpptraj_test_dir = amberhome + "/AmberTools/test/cpptraj/"
except (EnvironmentError, KeyError):
    amberhome = None
    has_amberhome = False
    cpptraj_test_dir = None

