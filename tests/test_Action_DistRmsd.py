import unittest; import pytraj as pt
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f1 = traj[1]
        act = adict['distrmsd']
        d0 = act('@CA', traj[:2], quick_get=True)
        print(d0[:])
        atommask = traj.top("@CA")
        print(f0.dist_rmsd(f1, traj.top('@CA')))
        assert f0.dist_rmsd(f1, traj.top('@CA')) == d0[1]
        print(f0.rmsd(f1, atommask))


if __name__ == "__main__":
    unittest.main()
