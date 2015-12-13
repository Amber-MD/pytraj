from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq



class Test(unittest.TestCase):

    def test_0(self):
        trajin = pt.datafiles.tc5b_trajin + """
        distance d0 @2 @3
        distance d1 @4 @7
        corr d0 d1 out test.out
        corr d0 out test2.out
        """

        # exclude DatasetTopology (1st)
        cout = pt.datafiles.load_cpptraj_output(trajin)[1:]

        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = pt.calc_distance(traj, ['@2 @3', '@4, @7'])

        pout = pt.xcorr(dslist[0],
                        dslist[1])
        # corr d0, d1
        aa_eq(pout, cout[2])

        # corr d0, d0
        pout = pt.xcorr(dslist[0],
                        dslist[0])
        aa_eq(pout, cout[3])

        # autocorr d0, d0
        pout = pt.acorr(dslist[0])
        aa_eq(pout, cout[3])


if __name__ == "__main__":
    unittest.main()
