from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestReplicateCell(unittest.TestCase):

    def test_vs_cpptraj(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        txt = '''
        parm ./data/tz2.ortho.parm7
        trajin ./data/tz2.ortho.nc
        replicatecell name test all
        '''

        t0 = pt.replicate_cell(traj, direction='all')
        state = pt.load_cpptraj_state(txt)
        state.run()
        saved_t0 = state.data[1]
        aa_eq(saved_t0.xyz, t0.xyz)

    def test_vs_list_tuple(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        traj0 = pt.replicate_cell(traj, direction='dir 001 dir 0-10')
        traj1 = pt.replicate_cell(traj, direction=('001', '0-10'))
        aa_eq(traj0.xyz, traj1.xyz)

        self.assertRaises(ValueError, lambda: pt.replicate_cell(traj, direction=traj[0]))


if __name__ == "__main__":
    unittest.main()
