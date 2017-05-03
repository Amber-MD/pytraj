from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.utils.context import capture_stdout

cm = '''
parm {}
trajin {}
dihedral phi2 :1@C :2@N :2@CA :2@C type phi
dihedral phi3 :2@C :3@N :3@CA :3@C type phi
analyze crankshaft angle phi2 phi3
'''.format(fn('tz2.parm7'), fn('tz2.nc'))


class TestCrank(unittest.TestCase):
    def test_crank(self):
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        dihedrals = pt.dihedral(
            traj, [':1@C :2@N :2@CA :2@C', ':2@C :3@N :3@CA :3@C'])
        state = pt.load_cpptraj_state(cm)
        with capture_stdout() as (out, _):
            state.run()
        data = pt.crank(dihedrals[0], dihedrals[1], mode='angle')
        assert out.read() == data


if __name__ == "__main__":
    unittest.main()
