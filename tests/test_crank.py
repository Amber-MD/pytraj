from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


cm = '''
parm data/tz2.parm7
trajin data/tz2.nc
dihedral phi2 :1@C :2@N :2@CA :2@C type phi
dihedral phi3 :2@C :3@N :3@CA :3@C type phi
analyze crankshaft angle phi2 phi3
'''


class TestCrank(unittest.TestCase):

    def test_crank(self):
        traj = pt.iterload("data/tz2.nc", "data/tz2.parm7")
        dihedrals = pt.dihedral(traj, [':1@C :2@N :2@CA :2@C', ':2@C :3@N :3@CA :3@C'])
        state = pt.load_cpptraj_state(cm)
        state.run()
        # TODO: assert please
        # cpptraj does not dump data to Dataset
        pt.crank(dihedrals[0], dihedrals[1], mode='angle')


if __name__ == "__main__":
    unittest.main()
