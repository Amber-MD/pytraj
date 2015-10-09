from __future__ import print_function
import pytraj as pt
import numpy as np
import unittest
from pytraj.hbonds import search_hbonds, search_nointramol_hbonds
from pytraj.testing import aa_eq
from pytraj.compat import izip as zip


# TODO: assert
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

if __name__ == "__main__":
    unittest.main()
