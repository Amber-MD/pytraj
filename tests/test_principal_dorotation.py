from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir


class TestPrincipalAxis(unittest.TestCase):

    def test_do_rotation(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        cm = '''
        principal * dorotation mass name pout
        createcrd myname
        '''
        state = pt.load_cpptraj_state(cm, traj)
        state.run()

        mut_traj_0 = traj[:]
        mut_traj_1 = traj[:]

        data = pt.principal_axes(mut_traj_0, mask='*', dorotation=True, mass=True)
        pt.align_principal_axis(mut_traj_1, mask='*', mass=True)

        aa_eq(data[0], state.data[1].values)
        aa_eq(data[1], state.data[2].values)
        aa_eq(state.data[-1].xyz, mut_traj_0.xyz)
        aa_eq(state.data[-1].xyz, mut_traj_1.xyz)

if __name__ == "__main__":
    unittest.main()
