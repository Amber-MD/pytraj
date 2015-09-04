from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    # TODO : add assert
    # what cpptraj dump in dslist?

    def test_0(self):
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        traj = mdio.iterload(fn, fn)

        dslist = DataSetList()
        act = adict['nastruct']
        act("", traj, dslist=dslist)
        act.print_output()

        d0 = dslist[0]

        for d0 in dslist:
            pass

        dsmall = dslist.get_dataset(dtype='float')
        ds_int = dslist.get_dataset(dtype='integer')

        import numpy as np
        dsmall = np.asarray(dsmall)

    def test_1(self):
        from pytraj.common_actions import nastruct
        fn = "./data/Test_NAstruct/adh026.3.pdb"
        traj = mdio.iterload(fn, fn)
        dslist = nastruct(traj, dtype='dataset')
        dsize = dslist.size

        # dummy loops
        for i in range(100):
            dslist = nastruct(traj, dtype='dataset')
            assert dslist.size == dsize


if __name__ == "__main__":
    unittest.main()
