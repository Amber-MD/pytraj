from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

'''Aim: check segmentation fault for Topology that does not have bond info'''


class Test(unittest.TestCase):

    @test_if_having("mdtraj")
    def test_0(self):
        import mdtraj as md
        from mdtraj.testing import get_fn
        from pytraj import common_actions as pyca
        t = md.load(get_fn("frame0.gro"))
        traj = io.load_mdtraj(t, autoconvert=False)

        # FIXME: segmentation fault: no radii info, calc_molsurf
        # traj.calc_molsurf()

        traj.calc_dssp()
        traj.rmsd(mode='cpptraj')
        traj.rmsd(0)
        md.rmsd(t, t, 0)


if __name__ == "__main__":
    unittest.main()
