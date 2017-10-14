import unittest

import pytraj as pt
from pytraj.testing import aa_eq

from utils import fn


class Test(unittest.TestCase):
    def test_0(self):
        mask = "@CA"
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        traj.top
        atm = traj.top(mask)
        n_selected_atoms = atm.n_atoms
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        aa_eq(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        aa_eq(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        assert (newtraj.xyz.shape == (traj.n_frames, n_selected_atoms, 3))

        # check if there is segmentation fault
        i = 0
        for frame in traj:
            i += 1

    def test_1(self):
        # why Trajectory is here? because I am lazy to move
        mask = "@CA"
        # creat Trajectory ( [:] )
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))[:]
        traj.top
        atm = traj.top(mask)
        newtraj = traj[atm]
        newtraj2 = traj[mask + ' :frame']
        aa_eq(newtraj2.xyz.flatten(), newtraj.xyz.flatten())
        aa_eq(newtraj2.xyz.flatten(), newtraj.xyz.flatten())


if __name__ == "__main__":
    unittest.main()
