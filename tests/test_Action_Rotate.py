import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        act = adict['rotate']
        f0 = traj[0]
        f1 = traj[1]
        f0cp = f0.copy()
        f0saved = f0.copy()
        f0cp2 = f0.copy()
        f1cp = f1.copy()
        print (f0cp.rmsd(f0))
        print (f1cp.rmsd(f1))
        # rotate all CA atoms (I know it's very dummy test)
        act("@CA x 60 y 120 z 50", f0, traj.top)
        act.do_action(f1)
        print (f0cp.rmsd(f0))
        print (f1cp.rmsd(f1))

        f0cp2.rotate(60, 120, 50, traj.top('@CA'))
        print ("f0cp2.rmsd(f0saved) = ", f0cp2.rmsd(f0saved))

        # do more
        # create mutable FrameArray
        farray = traj[:]
        # perform action on farray
        act.do_action(farray)

        for f1, f2 in zip(farray, traj):
            print (f1.rmsd(f2))

    def test_1(self):
        print ("test 1: assert")
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0cp = f0.copy()
        act = adict['rotate']
        f0.rotate(60, 120, 50, traj.top('@CA'))
        act("@CA x 60 y 120 z 50", f0cp, traj.top)
        fsaved = mdio.load("./data/rotated_frame0.x60y120z50.Tc5b.r", "./data/Tc5b.top")[0]
        print (f0.rmsd(fsaved))
        print (f0cp.rmsd(fsaved))
        assert f0.rmsd(fsaved) < 1E-3
        assert f0cp.rmsd(fsaved) < 1E-3

        # do_rotation
        from pytraj.common_actions import do_rotation
        # make mutable FrameArray object from TrajReadOnly object
        farray = traj[:]
        # perform rotation
        do_rotation("@CA x 60 y 120 z 50", farray)
        # assert
        assert farray[0].rmsd(fsaved) < 1E-3

if __name__ == "__main__":
    unittest.main()
