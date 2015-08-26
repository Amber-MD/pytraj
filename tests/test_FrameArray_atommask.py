import unittest; import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        print("test Trajectory")
        mask = "@CA"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        top = traj.top
        frame0 = traj[0].copy()
        print('mask @CA for frame')
        print(frame0[top('@CA')])
        print('-------------------')
        f1 = frame0.get_subframe("@CA", traj.top)
        print(f1.size)
        print(f1[:2])
        frame0.strip_atoms("!@CA", traj.top)
        print(frame0[:2])
        farray0 = traj['@CA :frame']

        assert_almost_equal(farray0[0].coords, frame0.coords)

        _farray = Trajectory()
        _farray.top = traj.top._modify_state_by_mask(traj.top(mask))
        print(top('@CA').n_atoms)
        for i, frame in enumerate(traj):
            print(frame[top('@CA')])
            _frame = frame.get_subframe(mask, traj.top)
            _farray.append(_frame)
        print(_farray[0])
        print(_farray[0, :2])

    def test_1(self):
        print("test TrajectoryIterator")
        mask = "@CA"
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj)
        top = traj.top
        frame0 = traj[0].copy()
        print('mask @CA for frame')
        print(frame0[top('@CA')])
        print('-------------------')
        f1 = frame0.get_subframe("@CA", traj.top)
        print(f1.size)
        print(f1[:2])
        frame0.strip_atoms("!@CA", traj.top)
        print(frame0[:2])
        farray0 = traj['@CA :frame']

        assert_almost_equal(farray0[0].coords, frame0.coords)

        _farray = Trajectory()
        _farray.top = traj.top._modify_state_by_mask(traj.top(mask))
        print(top('@CA').n_atoms)
        for i, frame in enumerate(traj):
            print(frame[top('@CA')])
            _frame = frame.get_subframe(mask, traj.top)
            _farray.append(_frame)
        print(_farray[0])
        print(_farray[0, :2])

    def test_2(self):
        print("test TrajectoryIterator")
        mask = "@CA"
        # TrajectoryIterator
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # Trajectory
        farray = traj[:]

        frame0 = traj[0].copy()
        frame0.strip_atoms("!@CA", traj.top)

        frame1 = traj['@CA :frame'][0]
        frame2 = farray['@CA :frame'][0]

        assert frame0.rmsd(frame1) < 1E-3
        assert frame2.rmsd(frame1) < 1E-3


if __name__ == "__main__":
    unittest.main()
