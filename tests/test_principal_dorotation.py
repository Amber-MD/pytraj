from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir
import pytraj.all_actions as pyca


class TestPrincipalAxis(unittest.TestCase):

    def test_do_rotation(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0]
        f0 = traj[0].copy()
        pyca.align_principal_axis(frame, "@CA", top=traj.top)

        # assert
        saved_frame = pt.iterload("./data/Tc5b.principal_dorotation.rst7",
                                  traj.top)[0]
        assert (frame.rmsd_nofit(saved_frame)) < 0.15

        cm = 'principal * dorotation mass name pout'
        state = pt.load_cpptraj_state(cm, traj)
        state.run()
        data = pt.principal_axes(traj, mask='*', dorotation=True, mass=True)
        aa_eq(data[0], state.data[1].values)
        aa_eq(data[1], state.data[2].values)


if __name__ == "__main__":
    unittest.main()
