from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0_0, f0_1 = traj[0].copy(), traj[0].copy()

        f1_0 = traj[1].copy()
        print("before tran_rot_tran: f0_0 and f0_1 must have the same coords")
        print(f0_0[0])
        print(f0_1[0])
        rmsd, mat, v1, v2 = f0_0.rmsd(f1_0, get_mvv=True)
        # f0_0 coords were updated  too
        print(rmsd, mat, v1, v2)

        print("after trans_rot_trans")
        f0_1.trans_rot_trans(v1, mat, v2)
        print(f0_0[0])
        print(f0_1[0])
        print(
            "f0_0 and f0_1 have different coords. Am I missing anything here?")

    def test_1(self):
        trajr = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = trajr[0]
        f1 = trajr[1]

        arr0 = f0[0]
        arr1 = f1[0]

        f0.rmsd(f1)

        # make sure that calling `rmsd` does not change coords of frame
        # TODO: add keyword `update_frame=True` for updating coords
        assert_almost_equal(arr0, f0[0])
        assert_almost_equal(arr1, f1[0])

if __name__ == "__main__":
    unittest.main()
