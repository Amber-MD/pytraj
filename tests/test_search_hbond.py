from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
from pytraj.externals.six import iteritems as items
from pytraj.compat import izip as zip

class Test(unittest.TestCase):
    #@no_test
    def test_0(self):
        from pytraj.hbonds import search_hbonds, search_nointramol_hbonds
        traj = mdio.iterload("./data/DPDP.nc", "./data/DPDP.parm7")
        print ('n_frames = %s' % traj.n_frames)
        dslist = search_hbonds(traj)
        print (dslist.keys())
        for key in dslist.keys():
            if 'UU' not in key:
                assert dslist[key].tolist().__len__() == traj.n_frames
        mydict = dslist.to_dict()
        mydict_np = dslist.to_dict(use_numpy=True)
        assert len(mydict.keys()) == dslist.size
        assert len(mydict_np.keys()) == dslist.size
        import numpy as np

        for key in mydict.keys():
            mydict[key] = np.asarray(mydict[key])
            aa_eq(mydict[key], mydict_np[key])

        dslist_b = search_nointramol_hbonds(traj)
        print (dslist_b.size, dslist_b.keys())

    def test_1(self):
        print ("test get indices from legends")
        import numpy as np
        from pytraj.hbonds import search_hbonds
        from pytraj.misc import from_legends_to_indices
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        import pytraj.common_actions as pyca
        ds = pyca.search_hbonds(traj, dtype='dataset')
        print (ds.size)
        print (ds.to_dict())
        d0 = ds.groupby("@")
        legends = d0.keys()
        print (legends)
        indices = np.asarray(from_legends_to_indices(legends, traj.top))
        print (indices)

        print (traj.n_frames)
        shape = [traj.n_frames, indices.shape[0]]
        arr0 = np.empty(shape)
        print (arr0.shape)
        for i, frame in enumerate(traj):
            arr0[i] = frame.calc_distance(indices)

        d = dict(zip(legends, zip(*arr0)))
        print (d)
        arr1 = pyca.calc_distance(traj, indices, n_frames=traj.n_frames)
        assert_almost_equal(arr0.flatten(), arr1.flatten())

    def test_2(self):
        # test memory error
        traj = mdio.load("./data/Tc5b.crd", "./data/Tc5b.top")
        dslist0 = traj.search_hbonds()
        expected_n_hbonds = 6
        assert dslist0.groupby("UU").values[0] == expected_n_hbonds
        assert traj.search_hbonds().groupby("UU").values[0] == expected_n_hbonds

if __name__ == "__main__":
    unittest.main()
