from __future__ import absolute_import, print_function
#from pytraj.base import *
#from pytraj import io as mdio
#from pytraj.utils.check_and_assert import assert_almost_equal
from .data_sample.load_sample_data import load_sample_data

def run_tests():
    load_sample_data()
    print ("OK")
