from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.analysis_dict import AnalysisDict
analdict = AnalysisDict()
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):

        print (analdict.keys())

if __name__ == "__main__":
    unittest.main()
