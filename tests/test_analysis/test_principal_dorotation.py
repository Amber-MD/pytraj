from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn

from pytraj.testing import aa_eq


class TestPrincipalAxis(unittest.TestCase):
    def test_principal_axes_and_lign_principal_axis(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        cm = '''
        principal * dorotation mass name pout
        createcrd myname
        '''
        state = pt.load_cpptraj_state(cm, traj)
        state.run()

        mut_traj_0 = traj[:]
        mut_traj_1 = traj[:]

        data = pt.principal_axes(
            mut_traj_0, mask='*', dorotation=True, mass=True)
        pt.align_principal_axis(mut_traj_1, mask='*', mass=True)

        aa_eq(data[0], state.data[1].values)
        aa_eq(data[1], state.data[2].values)
        aa_eq(state.data[-1].xyz, mut_traj_0.xyz)
        aa_eq(state.data[-1].xyz, mut_traj_1.xyz)


if __name__ == "__main__":
    unittest.main()
