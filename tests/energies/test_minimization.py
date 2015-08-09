from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca
from pytraj.utils import has_
from pytraj.testing import amberhome


class Test(unittest.TestCase):

    @test_if_path_exists(amberhome)
    def test_0(self):
        print (amberhome)
        from pytraj.amber_wrap import minimize

        traj = pt.iterload("./data/Ala3/Ala3.crd",
                  "./data/Ala3/Ala3.top")
        t0 = traj[:1]

        if has_("sander"):
            print ('egb: ', pt.energy_decomposition(t0, igb=8, verbose=False)['gb'])

        print (pt.rmsd(traj, ref=t0[0], mask='!@H='))
        minimize(t0)
        print (pt.rmsd(traj, ref=t0[0], mask='!@H='))

        self.assertRaises(ValueError, lambda: minimize(traj))

        if has_("sander"):
            print ('egb: ', pt.energy_decomposition(t0, igb=8, verbose=False)['gb'])

        # load saved file
        saved_coords = pt.load("./data/Ala3/min/min.r", traj.top).xyz
        aa_eq(t0.xyz, saved_coords, decimal=5)


if __name__ == "__main__":
    unittest.main()
