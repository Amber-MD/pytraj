import unittest
from pytraj.six_2 import izip
import os
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.common_actions import do_rotation

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
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
        # create mutable Trajectory
        farray = traj[:]
        # perform action on farray
        act.do_action(farray)

        for f1, f2 in izip(farray, traj):
            print (f1.rmsd(f2))

    def test_1(self):
        print ("test 1: assert")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        f0 = traj[0]
        f0cp = f0.copy()
        act = adict['rotate']
        f0.rotate(60, 120, 50, traj.top('@CA'))
        act("@CA x 60 y 120 z 50", f0cp, traj.top)
        fsaved = mdio.iterload("./data/rotated_frame0.x60y120z50.Tc5b.r", "./data/Tc5b.top")[0]
        print (f0.rmsd(fsaved))
        print (f0cp.rmsd(fsaved))
        assert f0.rmsd(fsaved) < 1E-3
        assert f0cp.rmsd(fsaved) < 1E-3

        # do_rotation
        from pytraj.common_actions import do_rotation
        # make mutable Trajectory object from TrajectoryIterator object
        farray = traj[:]
        farray_cp1 = farray.copy()
        farray_cp2 = farray.copy()
        # perform rotation
        do_rotation(farray, "@CA x 60 y 120 z 50")
        # assert
        assert farray[0].rmsd(fsaved) < 1E-3

        # test timing
        from pytraj.utils.Timer import Timer
        with Timer() as t1:
            do_rotation(farray_cp1, "@CA x 60 y 120 z 50")

        with Timer() as t2:
            for frame in farray_cp2:
                frame.rotate(60, 120, 50, traj.top('@CA'))

        # make sure get the same result
        for f1, f2 in izip(farray_cp1, farray_cp2):
            assert f1.rmsd(f2) < 1E-3

        print ("time t1 = %s (s), t2 = %s (s) " %(t1.time_gap, t2.time_gap))

    def test_2(self):
        if os.path.exists("./data/NuG2/test.x.000"):
            # do local test (not in travis) due to big file
            print ("test much bigger file")
            # test timing
            from pytraj.utils.Timer import Timer

            traj = mdio.iterload("./data/NuG2/test.x.000", "./data/NuG2/NuG2.top")
            farray = traj[:]
            farray_cp1 = farray.copy()
            farray_cp2 = farray.copy()

            with Timer() as t1:
                do_rotation("@CA x 60 y 120 z 50", farray_cp1)

            with Timer() as t2:
                for frame in farray_cp2:
                    frame.rotate(60, 120, 50, traj.top('@CA'))

            # make sure get the same result
            for f1, f2 in izip(farray_cp1, farray_cp2):
                assert f1.rmsd(f2) < 1E-3

            tgap1 = t1.time_gap
            tgap2 = t2.time_gap
            print ("time for do_rotation: t1 = %s (s), frame.rotate: t2 = %s (s) " 
                  %(t1.time_gap, t2.time_gap))
            print ("speed up using do_rotation vs frame.rotate = %s" % (tgap2/tgap1))

if __name__ == "__main__":
    unittest.main()
