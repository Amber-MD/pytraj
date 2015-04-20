from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
from pytraj.externals.six import iteritems as items

class Test(unittest.TestCase):
    @no_test
    def test_0(self):
        # TODO : need to check with DRR about the result
        from pytraj.hbonds import search_hbonds
        traj = mdio.load("./data/DPDP.nc", "./data/DPDP.parm7")
        print ('n_frames = %s' % traj.n_frames)
        dslist = search_hbonds(traj)
        print (dslist.keys())
        for key in dslist.keys():
            if 'UU' not in key:
                print (key, dslist[key].tolist())
        mydict = dslist.to_dict()
        assert len(mydict.keys()) == dslist.size

    def test_1(self):
        # TODO : need to check with DRR about the result
        from pytraj.hbonds import search_hbonds
        traj = mdio.load("./data/Tc5b.crd", "./data/Tc5b.top")[:1]
        import pytraj.common_actions as pyca
        ds = pyca.search_hbonds(traj)
        print (ds.size)
        print (ds.to_dict())

if __name__ == "__main__":
    unittest.main()
