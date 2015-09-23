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

        f0.rmsd(f1)
        assert_almost_equal(arr0, f0[0])
        assert_almost_equal(arr1, f1[0])

        f1.rmsfit(f0)

        # expect reference `f0` coords are not changed
        assert_almost_equal(arr0, f0[0])

        trajsaved = mdio.iterload(
            "./data/fit_to_1stframe.Tc5b.x", "./data/Tc5b.top")
        f1saved = trajsaved[1]

        # make sure we reproduce cpptraj output
        assert_almost_equal(f1.coords, f1saved.coords, decimal=3)

        farray = traj[:]
        farray.rmsfit(traj[0])
        assert_almost_equal(farray[1].coords, f1saved.coords, decimal=3)

        farray.rmsfit('first')


if __name__ == "__main__":
    unittest.main()
