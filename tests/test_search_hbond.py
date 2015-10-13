from __future__ import print_function
import pytraj as pt
import numpy as np
import unittest
from pytraj.hbonds import search_hbonds, search_hbonds_nointramol
from pytraj.testing import aa_eq
from pytraj.compat import izip as zip


class TestSearchHbonds(unittest.TestCase):
    def test_hbonds(self):
        traj = pt.iterload("./data/DPDP.nc", "./data/DPDP.parm7")
        dslist = search_hbonds(traj, dtype='dataset')
        for key in dslist.keys():
            if 'UU' not in key:
                assert dslist[key].tolist().__len__() == traj.n_frames
        mydict = dslist.to_dict()
        mydict_np = dslist.to_dict(use_numpy=True)
        assert len(mydict.keys()) == dslist.size
        assert len(mydict_np.keys()) == dslist.size

        for key in mydict.keys():
            mydict[key] = np.asarray(mydict[key])
            aa_eq(mydict[key], mydict_np[key])

        dslist_b = search_nointramol_hbonds(traj)

    def test_hbonds_with_image(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        hbonds_0 = pt.search_hbonds(traj(autoimage=True))
        hbonds_1 = pt.search_hbonds(traj, image=True)
        aa_eq(hbonds_0.values, hbonds_1.values)

    def test_hbonds_from_pdb(self):
        traj = pt.load('data/1L2Y.pdb')
        hb = pt.search_hbonds(traj)

        state = pt.load_cpptraj_state('''
        parm data/1L2Y.pdb
        trajin data/1L2Y.pdb
        hbond series
        ''')
        state.run()

        for data_p, data_cpp in zip(hb.data, state.data[1:]):
            assert len(data_p) == traj.n_frames == 38, 'size of dataset must be 38'
            aa_eq(data_p, data_cpp)


if __name__ == "__main__":
    unittest.main()
