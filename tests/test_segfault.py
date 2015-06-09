from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.decorators import no_test, test_if_having
import pytraj.common_actions as pyca

"""
try not to get segmentation fault error (due to whatever freaking reason)
"""
traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

class Test(unittest.TestCase):
    def test_0(self):
        print ("iter")
        from pytraj._shared_methods import _frame_iter_master
        it = _frame_iter_master(traj)

        for idx, frame in enumerate(it):
            pass
            #print ("segmentation fault if uncommenting #traj")
            # Status: don't need to fix since "it" and "traj" share the same iterator
            #traj[idx]

        fa = traj[:]
        for idx, frame in enumerate(fa):
            fa[idx]

    def test_1(self):
        print ("calling search_hbonds several times")
        import pytraj.common_actions as pyca
        pyca.search_hbonds(traj)
        pyca.search_hbonds(traj, 'series')
        pyca.search_hbonds(traj, 'series, nointramol')

    def test_2(self):
        print ("DataSetList lifetime")
        #d = pyca.search_hbonds(traj)
        d = pyca.search_hbonds(traj).groupby("SER")
        d2 = pyca.search_hbonds(traj).groupby("SER").to_ndarray()
        print (d.size)
        print (d.keys())
        print (d)
        print (d2)

    def test_3_vdw_radii_topology(self):
        top = io.load_pdb("./data/tz2.pdb").top
        # should raise ValueError since pdb does not have vdw info
        self.assertRaises(ValueError, lambda: top.vdw_radii())

    def test_4_trajiter(self):
        traj = io.load_sample_data("tz2")
        from pytraj.compat import zip

        for idx, (f0, f1) in enumerate(zip(traj, traj)):
            f0.rmsd(f1)
        print (idx)
        #assert idx == traj.n_frames

    def test_indexing_nonrefernce_DSL(self):
        from pytraj import dihedral_analysis as da
        from pytraj.hbonds import search_hbonds

        # segmentation fault
        print (da.calc_phi(traj)[0])
        print (search_hbonds(traj)[0])
        # new DSL
        print (search_hbonds(traj)[:][:][:][:][0])
        print (search_hbonds(traj)[:3][2][:][:][0])
        d0_dummy = search_hbonds(traj)[:][:][:][:][0]
        d0 = search_hbonds(traj)[0]
        aa_eq(d0_dummy.to_ndarray(), d0.to_ndarray())
        # groupby
        print (search_hbonds(traj)[0:4].groupby("").groupby(""))

        dslist = da.calc_phi(traj)
        x = dslist[0]
        print (x)
        print (da.calc_phi(traj, dtype='dict'))


if __name__ == "__main__":
    unittest.main()
