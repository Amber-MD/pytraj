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

        print(f0[0], f1[0])
        f0.rmsd(f1)
        print(f0[0], f1[0])

        print("test_frame_fit: after fit_to")
        f1.rmsfit(f0)
        print(f0[0], f1[0])

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
        print(rmsd_1, rmsd_1_nofit)
        assert rmsd_1 - rmsd_1_nofit < 1E-3

        print("test farray.fit_to")
        farray.rmsfit(f0)
        print(farray[1].rmsd_nofit(f0))
        print(farray[1].rmsd(f0))
        assert rmsd_1 - farray[1].rmsd_nofit(f0) < 1E-3
        print(farray.rmsfit)

    def test_1(self):
        print("compare to cpptraj")

        # load frames to immutable traj
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(traj[0, 0])
        trajsaved = mdio.iterload(
            "./data/fit_to_1stframe.Tc5b.x", "./data/Tc5b.top")

        print("test trajsaved coords")
        for _f1 in trajsaved:
            print(_f1[0])
        print("END test trajsaved coords")

        f0saved = traj[0].copy()
        first = traj[0].copy()

        # make mutable traj
        farray = traj[:]

        assert_almost_equal(farray[0].coords, first.coords)
        print(farray)
        print("before fitting")
        print(farray[0, 0])
        farray.rmsfit(first, "*")
        print("after fitting")
        print(farray[0, 0])

        print("make sure to reproduce cpptraj output")
        for i, _f0 in enumerate(farray):
            _f1 = trajsaved[i]
            print(_f0[0], _f1[0])
            assert_almost_equal(_f0.coords, _f1.coords)
        print("END")

if __name__ == "__main__":
    unittest.main()
