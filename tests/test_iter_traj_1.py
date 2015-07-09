import unittest
import numpy as np
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        for i, f0 in enumerate(traj):
            for j, x in enumerate(f0.coords):
                if np.abs(x - 5.707) < 1E-3:
                    print("catch index index %s in %s-frame" % (j, i))

        print(traj[0, 0])
        print(traj[-1, 0])

        print(traj)
        for frame in traj.frame_iter():
            print(frame[0])

        for frame0 in farray.frame_iter():
            print(frame0[0])

        i = 0
        print('test frame_iter with stride')
        print(farray.n_frames)
        for frame0 in farray.frame_iter(start=0, stride=1):
            print(frame0)
            i += 1

        print(i)
        print(frame0[:10])
        print(traj[-1, :10])
        assert_almost_equal(traj[-1].coords, frame0.coords)

        for frame0 in farray.frame_iter(start=2, stride=4, stop=8):
            pass
        assert_almost_equal(traj[6].coords, frame0.coords)

        for frame0 in farray.frame_iter(start=2, stride=2):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        arr0 = traj[6][0]
        print("traj[6][0]")
        print(arr0)
        for frame0 in traj.frame_iter(start=2, stride=4, stop=8):
            pass
        print("traj[6][0]")
        print(traj[6][0])

        for frame0 in traj.frame_iter(start=2, stride=4, stop=8):
            pass
        print("traj[6][0]")
        print(traj[6][0])
        print(traj[5][0])
        print('frame0', frame[0])
        assert_almost_equal(traj[6].coords, frame0.coords)

        for frame0 in traj.frame_iter(start=2, stride=2):
            pass
        assert_almost_equal(traj[8].coords, frame0.coords)

        count = 0
        for frame0 in traj.frame_iter(start=2):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[-1].coords, frame0.coords)

        count = 0
        for frame0 in traj.frame_iter(start=2, stop=7):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[6].coords, frame0.coords)

        for frame0 in traj.frame_iter():
            pass
        assert_almost_equal(traj[-1].coords, frame0.coords)

        for frame0 in farray.frame_iter():
            pass
        assert_almost_equal(traj[-1].coords, frame0.coords)

        for frame0 in traj():
            pass
        assert_almost_equal(traj[-1].coords, frame0.coords)

        for frame0 in farray():
            pass
        assert_almost_equal(farray[-1].coords, frame0.coords)

        count = 0
        for frame0 in traj(start=2, stop=7):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[6].coords, frame0.coords)

        count = 0
        for frame0 in farray(start=2, stop=7):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[6].coords, frame0.coords)

        count = 0
        for frame0 in farray(2, 7, 1):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[6].coords, frame0.coords)

        count = 0
        for frame0 in farray(2, 7, 2):
            count += 1
            pass
        print('count = ', count)
        assert_almost_equal(traj[6].coords, frame0.coords)

    def test_1(self):
        print("test frame_iter with mask")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        top2 = traj.top.copy()

        for frame in traj(mask='@CA'):
            pass
        top2.strip_atoms("!@CA")
        assert frame.n_atoms == top2.n_atoms

        for frame in traj():
            pass
        assert frame.n_atoms == traj[0].n_atoms

        for frame in traj[:](mask='@CA'):
            pass

        f0 = traj[-1]
        f0.strip_atoms('!@CA', traj.top)
        assert_almost_equal(f0.coords, frame.coords)
        print(frame.n_atoms)
        assert frame.n_atoms == top2.n_atoms

        for frame in traj[:]():
            pass
        print(frame.n_atoms)
        assert frame.n_atoms == traj[0].n_atoms

    def test_2(self):
        print("test frame_iter with mask")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        for frame in farray(mask='@CA'):
            pass
        f0 = traj[-1]
        f0.strip_atoms('!@CA', traj.top)
        print(frame.n_atoms)
        assert_almost_equal(f0.coords, frame.coords)

if __name__ == "__main__":
    unittest.main()
