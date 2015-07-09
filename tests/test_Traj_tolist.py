from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    @test_if_having("numpy")
    def test_0(self):
        from numpy.testing import assert_almost_equal
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        xyz = traj.tolist()
        assert len(xyz) == traj.n_frames
        assert len(xyz[0]) == traj.n_atoms
        assert_almost_equal(xyz, traj.xyz)

if __name__ == "__main__":
    unittest.main()
