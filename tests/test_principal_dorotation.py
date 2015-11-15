from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame = traj[0]
        f0 = traj[0].copy()
        pyca.align_principal_axis(frame, "@CA", top=traj.top)

        # assert
        saved_frame = mdio.iterload("./data/Tc5b.principal_dorotation.rst7",
                                    traj.top)[0]
        assert (frame.rmsd_nofit(saved_frame)) < 0.15


if __name__ == "__main__":
    unittest.main()
