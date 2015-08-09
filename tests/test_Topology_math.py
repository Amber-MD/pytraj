from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists


class Test(unittest.TestCase):
    def test_0(self):
        t0 = mdio.load_sample_data().top
        t1 = mdio.load_sample_data('tz2').top

        t2 = t0 + t1  # mimic ParmEd
        assert t2.n_atoms == t0.n_atoms + t1.n_atoms

        t0 += t1
        assert t0.n_atoms == t2.n_atoms
        mdio.write_parm("output/test_join_topologies.prmtop", t0)


if __name__ == "__main__":
    unittest.main()
