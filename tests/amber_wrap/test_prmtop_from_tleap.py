from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca

class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.testing import amberhome

        if amberhome:
            from pytraj.amber_wrap import prmtop_from_tleap
            t0 = prmtop_from_tleap('./data/tz2.pdb')
            t1 = pt.load_topology('./data/tz2.pdb')
            print (t0, t1)
            assert t0.n_atoms == t1.n_atoms
        else:
            print ("does not have AMBERHOME. skip")

if __name__ == "__main__":
    unittest.main()
