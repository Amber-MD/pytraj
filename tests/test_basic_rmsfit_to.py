import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_frame_fit(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f1 = traj[1]

        arr0 = list(f0[0])
        arr1 = list(f1[0])

        print("before calling rmsd method")
        print(f0[0], f1[0])
        f0.rmsd(f1)
        print("after calling rmsd method: exect coords does not change")
        print(f0[0], f1[0])
        assert_almost_equal(arr0, f0[0])
        assert_almost_equal(arr1, f1[0])

        print("test_frame_fit: after fit_to: exect coords will be updated")
        f1.rmsfit(f0)
        print(f0[0], f1[0])

        # expect reference `f0` coords are not changed
        assert_almost_equal(arr0, f0[0])

        trajsaved = mdio.iterload(
            "./data/fit_to_1stframe.Tc5b.x", "./data/Tc5b.top")
        f1saved = trajsaved[1]
        print("from cpptraj, f1saved")
        print(f1saved[0])

        # make sure we reproduce cpptraj output
        assert_almost_equal(f1.coords, f1saved.coords)

        print("test: mutable traj (farray)")
        farray = traj[:]
        print(farray)
        print(farray[1, 0])
        farray.rmsfit(traj[0])
        print(farray[1, 0])
        assert_almost_equal(farray[1].coords, f1saved.coords)

        farray.rmsfit('first')


if __name__ == "__main__":
    unittest.main()
