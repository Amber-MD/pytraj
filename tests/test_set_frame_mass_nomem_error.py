from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        t0 = traj[:]
        f0 = t0[0]
        aa_eq(f0.xyz, traj[0].xyz)
        f0.set_frame_mass(traj.top)

        # make sure that frame's coords does not 
        # change after setting mass
        #print(f0.xyz)
        aa_eq(f0.xyz, traj[0].xyz)


if __name__ == "__main__":
    unittest.main()
