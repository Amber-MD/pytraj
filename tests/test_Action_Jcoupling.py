import os
import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj import adict
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.core.DataFileList import DataFileList
from pytraj.datasets.DataSetList import DataSetList
from pytraj.testing import cpptraj_test_dir, aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        kfile = os.path.join(cpptraj_test_dir, "Test_Jcoupling", "Karplus.txt")
        traj = mdio.iterload("./data/tz2.nc", "./data/tz2.parm7")
        frame = traj[0]
        command = "kfile %s" % kfile
        dslist = DataSetList()
        dflist = DataFileList()
        d0 = adict['jcoupling'](command, frame, traj.top,
                                dslist=dslist,
                                dflist=dflist)

        # another way, and assert
        dslist = DataSetList()
        from pytraj.common_actions import calc_jcoupling
        d0 = calc_jcoupling(traj, command)
        d1 = calc_jcoupling(traj, kfile=kfile)
        aa_eq(d0.values, d1.values)


if __name__ == "__main__":
    unittest.main()
