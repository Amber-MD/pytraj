from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Testtleap_wrapper(unittest.TestCase):

    def test_tleap(self):
        from pytraj.testing import amberhome

        if amberhome and os.path.exists(amberhome + '/bin/tleap'):
            from pytraj.amber_wrapper import prmtop_from_tleap
            t0 = prmtop_from_tleap('./data/tz2.pdb')
            t1 = pt.load_topology('./data/tz2.pdb')
            print(t0, t1)
            assert t0.n_atoms == t1.n_atoms
        else:
            print("does not have AMBERHOME. skip")


if __name__ == "__main__":
    unittest.main()
