from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq

"""
try not to get segmentation fault error (due to whatever freaking reason)
"""
traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.shared_methods import iterframe_master
        it = iterframe_master(traj)

        for idx, frame in enumerate(it):
            pass
            # Status: don't need to fix since "it" and "traj" share the same iterator

        fa = traj[:]
        for idx, frame in enumerate(fa):
            fa[idx]

    def test_1(self):

        pt.search_hbonds(traj)
        pt.search_hbonds(traj, 'series')
        pt.search_hbonds(traj, 'series, nointramol')

    def test_4_trajiter(self):
        traj = pt.load_sample_data("tz2")
        from pytraj.compat import zip

        for idx, (f0, f1) in enumerate(zip(traj, traj)):
            f0.rmsd(f1)

    def test_indexing_nonrefernce_DSL(self):
        from pytraj import dihedral_analysis as da

        # segmentation fault
        # new DSL
        d0_dummy = pt.search_hbonds(traj)[:][:][:][:][0]
        d0 = pt.search_hbonds(traj)[0]
        # filter

        dslist = da.calc_phi(traj)
        x = dslist[0]


if __name__ == "__main__":
    unittest.main()
