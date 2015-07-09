from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.misc import get_atts
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(get_atts(traj))

        dslist = DataSetList()
        act = adict['radgyr']
        act("", traj, dslist=dslist)

        alist = get_atts(dslist)
        print([a for a in alist if a.startswith("get_")])

        d0 = ['get_dtypes', 'get_scalar_modes', 'get_scalar_types',
              'get_aspects', 'get_legends']

        for att_str in d0:
            print(getattr(dslist, att_str)())

        print(dslist[0][:])
        print(dslist[1][:])

        top = getattr(traj, 'top')
        print(top)


if __name__ == "__main__":
    unittest.main()
