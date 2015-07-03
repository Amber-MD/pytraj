from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        it0, it1 = traj.make_independent_iterators(2)

        for idx, frame in enumerate(it0):
            pass
        assert idx == traj.n_frames - 1

        for idx, frame in enumerate(it1):
            pass
        assert idx == traj.n_frames - 1

        from pytraj.compat import zip
        it0, it1 = traj.make_independent_iterators(2)
        for idx, (f0, f1) in enumerate(zip(it0, it1)):
            assert f0.rmsd(f1) < 1E-5
        assert idx == traj.n_frames - 1


if __name__ == "__main__":
    unittest.main()
