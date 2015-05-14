from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.utils import Timer, has_
from pytraj.externals.six.moves import range

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        f0 = traj[0]
        f0.strip_atoms(mask='!@CA', top=traj.top)
        atm0 = traj.top("!@CA")
        atm0.invert_mask()
        NATOM = atm0.n_atoms
        assert f0.n_atoms == NATOM
        print (f0)

        fa0 = traj[:]
        fa0.strip_atoms('!@CA')
        assert fa0[0].n_atoms == NATOM
        fa1 = traj[:]
        fa1._fast_strip_atoms('!@CA')
        assert fa1[0].n_atoms == NATOM

        aa_eq(fa0.xyz, fa1.xyz)

        #fa3 = traj[:]
        fa3 = traj._fast_slice(slice(2, 100, 10))
        print (traj)
        #print (fa3)
        for i in range(10):
            fa3 += fa3.copy()

        saved_fa3 =  fa3.copy()

        fa4 = fa3.copy()
        print (fa3, fa4)

        @Timer()
        def normal_strip():
            fa3.strip_atoms("!@H,C")

        @Timer()
        def fast_strip_atoms():
            fa4._fast_strip_atoms("!@H,C")

        fa5 = fa3.copy()
        @Timer()
        def fancy_indexing():
            # as fast as _fast_strip_atoms
            fa5["@H,C"]

        print ("normal_strip")
        normal_strip()
        print ("fast_strip_atoms")
        fast_strip_atoms()
        print ("fancy_indexing")
        fancy_indexing()
        _fa5 = fa5['@H,C']
        aa_eq(fa3.xyz, fa4.xyz)
        aa_eq(_fa5.xyz, fa4.xyz)
        atm = traj.top("!@H,C")
        atm.invert_mask()
        assert fa3.n_atoms == atm.n_atoms == fa4.n_atoms == _fa5.n_atoms

        if has_("mdtraj"):
            import mdtraj as md
            m_traj_orig = md.load(traj.filename, top="./data/tz2.ortho.parm7")
            m_traj = m_traj_orig[2:100:10]
            m_traj.xyz = m_traj.xyz * 10. # convert to angstrom

            for i in range(10):
                m_traj += m_traj
            assert m_traj.n_frames == fa3.n_frames

            aa_eq(saved_fa3.xyz, m_traj.xyz)

            indices = saved_fa3.top('@H,C').indices
            s = m_traj.atom_slice(indices)
            aa_eq(fa3.xyz, s.xyz)

            @Timer()
            def slice_mdtraj_precalculated_indices():
                m_traj.atom_slice(indices)
            print ("slice_mdtraj_precalculated_indices:")
            slice_mdtraj_precalculated_indices() # ~6 times slower than pytraj's fa['@H,C']

            # remove solvent
            fa5 = saved_fa3.copy()

            @Timer()
            def pytraj_strip_wat():
                fa5['!:WAT']

            @Timer()
            def mdtraj_remove_solvent():
                m_traj.remove_solvent()

            fa5new = fa5['!:WAT']
            m_new = m_traj.remove_solvent()
            aa_eq(fa5new.xyz, m_new.xyz)

            print ("pytraj_strip_wat")
            pytraj_strip_wat()
            print ("mdtraj_remove_solvent")
            mdtraj_remove_solvent()
        else:
            print ("does not have mdtraj, skip")

if __name__ == "__main__":
    unittest.main()
