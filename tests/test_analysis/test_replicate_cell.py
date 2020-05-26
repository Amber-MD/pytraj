from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import pytest


class TestReplicateCell(unittest.TestCase):
    def test_vs_cpptraj(self):
        traj = pt.iterload(
            fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'), frame_slice=[(0, 1)])
        txt = '''
        parm {}
        trajin {} 1 1
        replicatecell name test dir 001
        '''.format(fn('tz2.ortho.parm7'), fn('tz2.ortho.nc'))

        t0 = pt.replicate_cell(traj, direction='dir 001')
        state = pt.load_cpptraj_state(txt)
        state.run()
        saved_t0 = state.data[1]
        aa_eq(saved_t0.xyz, t0.xyz)

    def test_vs_list_tuple(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        traj0 = pt.replicate_cell(traj, direction='dir 001 dir 0-10')
        traj1 = pt.replicate_cell(traj, direction=('001', '0-10'))
        aa_eq(traj0.xyz, traj1.xyz)

        with pytest.raises(ValueError):
            pt.replicate_cell(traj, direction=traj[0])


if __name__ == "__main__":
    unittest.main()
