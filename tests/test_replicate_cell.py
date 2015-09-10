from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class TestReplicateCell(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
        txt = '''
        parm ./data/tz2.ortho.parm7
        trajin ./data/tz2.ortho.nc
        replicatecell name test all
        '''
        t0 = pt.replicate_cell(traj, direction='all')
        cpp_out = pt.datafiles.load_cpptraj_output(txt, dtype='state')
        saved_t0 = cpp_out.datasetlist[0]
        aa_eq(saved_t0.xyz, t0.xyz)


if __name__ == "__main__":
    unittest.main()
