from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.utils import has_
from pytraj.testing import amberhome

try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False

has_sander_exe = os.path.exists(amberhome + '/bin/sander')
print('has_sander_exe', has_sander_exe)


class TestMin(unittest.TestCase):

    @unittest.skipIf(not has_sander or not has_sander_exe, 'does not have libsander/sander. skip')
    def test_0(self):
        if os.path.exists(amberhome):
            from pytraj.amber_wrapper import minimize

            traj = pt.iterload("./data/Ala3/Ala3.crd", "./data/Ala3/Ala3.top")
            t0 = traj[:1]

            if has_("sander"):
                print('egb: ', pt.energy_decomposition(t0, igb=8)['gb'])

            minimize(t0)

            self.assertRaises(ValueError, lambda: minimize(traj))

            if has_("sander"):
                print('egb: ', pt.energy_decomposition(t0, igb=8)['gb'])

            # load saved file
            saved_coords = pt.load("./data/Ala3/min/min.r", traj.top).xyz
            aa_eq(t0.xyz, saved_coords, decimal=5)


if __name__ == "__main__":
    unittest.main()
