from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0_0 = traj[0].copy()
        f0_1 = traj[0].copy()
        f0_2 = traj[0].copy()
        f0_3 = traj[0].copy()
        f0_4 = traj[0].copy()
        f1_0 = traj[1].copy()
        f1_1 = traj[1].copy()

        print ("1st try: fit to f1_0")
        print ("before calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])
        f0_0.rmsd(f1_0)
        print ("after calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])

        print  ()
        print ("2nd try: fit to f1_0")
        print ("before calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])
        f0_0.rmsd(f1_0)
        print ("after calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])

        print  ()
        print ("3rd try: fit to f1_1")
        print ("before calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])
        f0_0.rmsd(f1_1)
        print ("after calling `rmsd`")
        print (f0_0[0])
        print (f1_0[0])
        print (f1_1[0])


        f1_2 = traj[1].copy()
        print ("before tran_rot_tran: f0_3 and f0_4 must have the same coords")
        print (f0_3[:2])
        print (f0_4[:2])
        rmsd, mat, v1, v2 = f0_3.rmsd(f1_2, get_mvv=True)
        # f0_3 coords were updated  too
        print ('f0_3 updated:', f0_3[:2])
        print (rmsd, mat, v1, v2)

        print ("after trans_rot_trans")
        f0_4.trans_rot_trans(v1, mat, v2)
        print (f0_3[:2])
        print (f0_4[:2])
        print ("f0_3 and f0_4 have different coords. Am I missing anything here?")

if __name__ == "__main__":
    unittest.main()
