from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.compat import zip

traj = pt.iterload(top="./data/Tc5b.top",
                   filename='data/md1_prod.Tc5b.x', )

class TestGeometry(unittest.TestCase):
    def testRadgyr(self):
        txt = '''
        parm ./data/Tc5b.top
        trajin ./data/md1_prod.Tc5b.x
        radgyr @CA nomax
        radgyr nomax
        radgyr !@H= nomax
        '''

        data = pt.datafiles.load_cpptraj_output(txt)
        for mask, out  in zip(['@CA', '', '!@H='], data):
            aa_eq(pt.radgyr(traj, mask), out)


if __name__ == "__main__":
    unittest.main()
