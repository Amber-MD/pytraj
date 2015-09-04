import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_1(self):
        mask = "@CA"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top = traj.top
        frame0 = traj[0].copy()
        f1 = frame0.get_subframe("@CA", traj.top)
        frame0.strip_atoms("!@CA", traj.top)
        farray0 = traj['@CA ']

        assert_almost_equal(farray0[0].coords, frame0.coords)

        _farray = Trajectory()
        _farray.top = traj.top._modify_state_by_mask(traj.top(mask))
        for i, frame in enumerate(traj):
            _frame = frame.get_subframe(mask, traj.top)
            _farray.append(_frame)

    def test_2(self):
        mask = "@CA"
        # TrajectoryIterator
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # Trajectory
        farray = traj[:]

        frame0 = traj[0].copy()
        frame0.strip_atoms("!@CA", traj.top)

        frame1 = traj['@CA'][0]
        frame2 = farray['@CA'][0]

        assert frame0.rmsd(frame1) < 1E-3
        assert frame2.rmsd(frame1) < 1E-3


if __name__ == "__main__":
    unittest.main()
