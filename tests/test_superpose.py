import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class TestBasic(unittest.TestCase):
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

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        f0 = traj[0]
        f0saved = f0.copy()
        f1 = traj[1]

        rmsd_0 = f0.rmsd(f1)
        rmsd_0_nofit = f0.rmsd_nofit(f1)
        assert rmsd_0 != rmsd_0_nofit

        # do fitting
        f1.rmsfit(f0)
        rmsd_1 = f1.rmsd(f0)
        rmsd_1_nofit = f1.rmsd_nofit(f0)

        # make sure that rmsd_nofit after do fitting is equal to rmsd (with
        # fit)
        assert rmsd_1 - rmsd_1_nofit < 1E-3

        farray.rmsfit(f0)
        assert rmsd_1 - farray[1].rmsd_nofit(f0) < 1E-3

    def test_1(self):

        # load frames to immutable traj
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajsaved = mdio.iterload(
            "./data/fit_to_1stframe.Tc5b.x", "./data/Tc5b.top")

        for _f1 in trajsaved:
            pass

        f0saved = traj[0].copy()
        first = traj[0].copy()

        # make mutable traj
        farray = traj[:]

        assert_almost_equal(farray[0].coords, first.coords)
        farray.rmsfit(first, "*")
        farray2 = traj[:]
        farray2.superpose(first, "*")

        for i, _f0 in enumerate(farray):
            _f1 = trajsaved[i]
            assert_almost_equal(_f0.coords, _f1.coords, decimal=3)

        for i, _f0 in enumerate(farray2):
            _f1 = trajsaved[i]
            assert_almost_equal(_f0.coords, _f1.coords, decimal=3)


if __name__ == "__main__":
    unittest.main()
