from __future__ import print_function
#import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test, test_if_having


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.compat import izip as zip
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj2 = mdio.iterload("./data/md.crazy_ext", "./data/Tc5b.top")

        for idx, (f0, f1) in enumerate(zip(traj, traj2)):
            assert_almost_equal(f0.coords, f1.coords)
        assert idx + 1 == traj.n_frames == traj2.n_frames

if __name__ == "__main__":
    unittest.main()
