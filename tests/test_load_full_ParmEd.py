from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having

class Test(unittest.TestCase):
    @test_if_having("parmed")
    def test_0(self):
        top_fn = "./data/Tc5b.top"
        p_top = mdio._load_parmed(top_fn)
        top = mdio.load_full_ParmEd(p_top)
        assert len(p_top.atoms) == top.n_atoms

        p_top = mdio._load_parmed("./data/tz2.pdb")
        print (p_top)
        top = mdio.load_full_ParmEd(p_top)
        assert top.filename == "mytmptop.pdb"
        assert len(p_top.atoms) == top.n_atoms


if __name__ == "__main__":
    unittest.main()
