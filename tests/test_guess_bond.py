from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    @test_if_having("mdtraj")
    def test_0(self):
        from mdtraj.testing import get_fn
        import mdtraj as md
        t = md.load(get_fn("frame0.gro"))
        traj = io.load_mdtraj(t)
        print(traj.calc_dssp(dtype='ndarray'))
        print(traj.calc_multidihedral("phi"))

if __name__ == "__main__":
    unittest.main()
