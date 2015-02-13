import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        print (traj[0, 0])
        print (traj[-1, 0])

        print (traj)
        for frame in traj.frame_iter():
            print (frame[0])

        for frame0 in farray.frame_iter():
            print (frame0[0])

        i = 0
        print ('test frame_iter with stride')
        print (farray.n_frames)
        for frame0 in farray.frame_iter(start=0, stride=0):
            print (frame0)
            i += 1

        print (i)
        print (frame0[:10])
        print (traj[-1, :10])
        assert_almost_equal(traj[-1].coords, frame0.coords)

        for frame0 in farray.frame_iter(start=2, stride=4, stop=8):
            pass
        assert_almost_equal(traj[6].coords, frame0.coords)

        for frame0 in farray.frame_iter(start=2, stride=2):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        for frame0 in traj.frame_iter(start=2, stride=4, stop=8):
            pass
        assert_almost_equal(traj[6].coords, frame0.coords)

        for frame0 in traj.frame_iter(start=2, stride=2):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        count = 0
        for frame0 in traj.frame_iter(start=2):
            count += 1
            pass
        print ('count = ', count)
        assert_almost_equal(traj[-1].coords, frame0.coords)

        count = 0
        for frame0 in traj.frame_iter(start=2, stop=7):
            count += 1
            pass
        print ('count = ', count)
        assert_almost_equal(traj[7].coords, frame0.coords)

        for frame0 in traj.frame_iter(indices=(1, 3, 4, 8)):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        for frame0 in farray.frame_iter(indices=(1, 3, 4, 5)):
            pass
        assert_almost_equal(traj[5].coords, frame0.coords)

        for frame0 in traj(indices=(1, 3, 4, 8)):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        for frame0 in farray(indices=(1, 3, 4, 5)):
            pass
        assert_almost_equal(traj[5].coords, frame0.coords)

        count = 0
        for frame0 in traj(start=2, stop=7):
            count += 1
            pass
        print ('count = ', count)
        assert_almost_equal(traj[7].coords, frame0.coords)

        count = 0
        for frame0 in farray(start=2, stop=7):
            count += 1
            pass
        print ('count = ', count)
        assert_almost_equal(traj[7].coords, frame0.coords)

if __name__ == "__main__":
    unittest.main()
