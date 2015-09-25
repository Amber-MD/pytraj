from __future__ import print_function
import unittest
import numpy as np
import pytraj as pt
from pytraj import Frame
from pytraj.utils import assert_almost_equal, Timer
from pytraj.utils import has_


class Test(unittest.TestCase):
    def test_xyz(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = Frame()
        frame.append_xyz(traj[0].xyz)
        assert_almost_equal(frame.coords, traj[0].coords)
        assert_almost_equal(frame.xyz.flatten(), traj[0].xyz.flatten())
        assert_almost_equal(np.array(frame.buffer1d), traj[0].xyz.flatten())

    @unittest.skipIf(not has_('mdtraj'), 'does not have mdtraj')
    def test_1(self):
        from mdtraj.formats import psf
        import mdtraj as md
        from mdtraj.testing import get_fn
        import numpy as np

        fname = get_fn('ala_ala_ala.pdb')
        m_traj = md.load(fname)
        #print(m_traj)
        f0 = Frame()
        f1 = f0.copy()
        f0.append_xyz(m_traj.xyz[0].astype(np.float64))
        farray = pt.load_mdtraj(m_traj, autoconvert=False, top=fname)
        f1 = farray[0]

        assert_almost_equal(f0.coords, f1.coords)


if __name__ == "__main__":
    unittest.main()
