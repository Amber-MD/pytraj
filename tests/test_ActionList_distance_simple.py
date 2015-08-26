from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]

        trajin = pt.datafiles.tc5b_trajin + """
        distance @CB @CA
        distance @CA @H
        """

        cout = pt.datafiles.load_cpptraj_output(trajin)
        #print("cout", cout)

        mask_list = ('@CB @CA', '@CA @H')
        dslist = pt.calc_distance(traj, mask_list)

        # compare to cpptraj output
        aa_eq(dslist.flatten(), cout.values.flatten())

        #print("@CB @CA", pt.calc_distance(traj, "@CB @CA"))

        dslist3_0 = pt.calc_distance(traj, mask_list[0])
        dslist3_1 = pt.calc_distance(traj, mask_list[1])

        aa_eq(dslist3_0, dslist[0])
        aa_eq(dslist3_1, dslist[1])


if __name__ == "__main__":
    unittest.main()
