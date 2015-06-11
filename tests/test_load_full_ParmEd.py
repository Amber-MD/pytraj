from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having("chemistry")
    def test_0(self):
        top_fn = "./data/Tc5b.top"
        p_top = mdio._load_parmed(top_fn)
        print (p_top)
        top = mdio.load_full_ParmEd(p_top)
        print (top)

if __name__ == "__main__":
    unittest.main()
