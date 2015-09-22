from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
import pytraj.common_actions as pyca
"""
try not to get segmentation fault error (due to whatever freaking reason)
"""
traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj._shared_methods import iterframe_master
        it = iterframe_master(traj)

        for idx, frame in enumerate(it):
            pass
            # Status: don't need to fix since "it" and "traj" share the same iterator
            # traj[idx]

        fa = traj[:]
        for idx, frame in enumerate(fa):
            fa[idx]

    def test_1(self):
        import pytraj.common_actions as pyca
        pyca.search_hbonds(traj)
        pyca.search_hbonds(traj, 'series')
        pyca.search_hbonds(traj, 'series, nointramol')

    def test_2(self):
        #d = pyca.search_hbonds(traj)
        d = pyca.search_hbonds(traj).filter("SER")
        d2 = pyca.search_hbonds(traj).filter("SER").to_ndarray()

    def test_3_vdw_radii_topology(self):
        top = io.load_pdb("./data/tz2.pdb").top
        # should raise ValueError since pdb does not have vdw info
        self.assertRaises(ValueError, lambda: top.vdw_radii())

    def test_4_trajiter(self):
        traj = io.load_sample_data("tz2")
        from pytraj.compat import zip

        for idx, (f0, f1) in enumerate(zip(traj, traj)):
            f0.rmsd(f1)
        #assert idx == traj.n_frames

    def test_indexing_nonrefernce_DSL(self):
        from pytraj import dihedral_analysis as da
        from pytraj.hbonds import search_hbonds

        # segmentation fault
        # new DSL
        d0_dummy = search_hbonds(traj)[:][:][:][:][0]
        d0 = search_hbonds(traj)[0]
        aa_eq(d0_dummy.to_ndarray(), d0.to_ndarray())
        # filter

        dslist = da.calc_phi(traj)
        x = dslist[0]


if __name__ == "__main__":
    unittest.main()
