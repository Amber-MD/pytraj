from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestDihedral(unittest.TestCase):

    def test_dihedral(self):
        import numpy as np
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        fa = traj[:]
        mask = ':2@CA :14@CA :15@CA :16@CA'
        txt = '''
        parm ./data/Tc5b.top
        trajin ./data/Tc5b.x
        dihedral %s
        ''' % mask
        d0 = pt.dihedral(traj, mask, dtype='dataset').to_ndarray()
        d1 = pt.dihedral(traj, mask)
        d2 = pt.calc_dihedral(fa, mask)
        state = pt.load_cpptraj_state(txt)
        state.run()
        dcpp = state.data[1:].values

        aa_eq(d0, d1)
        aa_eq(d0, d2)
        aa_eq(d0, dcpp)

        Nsize = 10
        np.random.seed(1)
        arr = np.random.randint(0, 300, size=Nsize * 4).reshape(Nsize, 4)
        d3 = pt.calc_dihedral(fa, arr)
        d4 = pt.dihedral(traj, arr)
        d5 = pt.dihedral(traj, arr)
        d6 = pt.dihedral(fa, arr)
        d7 = pt.dihedral([fa, traj], arr, n_frames=2 * fa.n_frames)
        aa_eq(d3, d4)
        aa_eq(d3, d5)
        aa_eq(d3, d6)
        aa_eq(d3.T, d7.T[:fa.n_frames])
        aa_eq(d3.T, d7.T[fa.n_frames:])

        d8 = pt.dihedral(traj, mask, dtype='dataset')
        d9 = pt.tools.dict_to_ndarray(pt.dihedral(traj, mask, dtype='dict'))
        aa_eq(d0, d8.values)
        aa_eq(d0, d9)

        # raise
        self.assertRaises(ValueError, lambda: pt.dihedrals(traj, [[0, 3, 2]]))


if __name__ == "__main__":
    unittest.main()
