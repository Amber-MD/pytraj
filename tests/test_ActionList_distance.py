from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
mask_list = ('@CB @CA', '@CA @H')


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.core.ActionList import ActionList
        from pytraj.datasets.DatasetList import DatasetList
        from pytraj.actions import CpptrajActions as CA

        dslist = DatasetList()
        actlist = ActionList()

        for mask in mask_list:
            actlist.add_action(
                CA.Action_Distance(), mask, traj.top,
                dslist=dslist)
        actlist.do_actions(traj)

        dslist2 = pyca.calc_distance(traj, mask_list)
        aa_eq(dslist.values, dslist2)

    def test_1(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        mask_list = ('@CB @CA', '@CA @H')
        dslist = pt.calc_distance(traj, mask_list)
        dslist3_0 = pt.calc_distance(traj, mask_list[0])
        dslist3_1 = pt.calc_distance(traj, mask_list[1])

        #print(dslist)
        #print(dslist3_0)
        #print(dslist3_1)

        aa_eq(dslist3_0, dslist[0])
        aa_eq(dslist3_1, dslist[1])


if __name__ == "__main__":
    unittest.main()
