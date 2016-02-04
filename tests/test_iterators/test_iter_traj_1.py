import unittest
import pytraj as pt
import numpy as np
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.compat import zip
from pytraj.testing import aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]

        for i, f0 in enumerate(traj):
            for j, x in enumerate(f0.xyz):
                if np.all(np.abs(x - 5.707) < 1E-3):
                    pass

        for frame in traj.iterframe():
            pass

        for frame0 in farray.iterframe():
            pass

        i = 0
        for frame0 in farray.iterframe(start=0, step=1):
            i += 1

        assert_almost_equal(traj[-1].xyz, frame0.xyz)

        start, stop, step = 2, 8, 4
        indices = list(range(start, stop, step))

        for idx, frame0, f in zip(indices, farray.iterframe(start, stop, step),
                                  traj[indices]):
            aa_eq(frame0.xyz, f.xyz)
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        for frame0 in farray.iterframe(start=2, step=2):
            pass
        assert_almost_equal(traj[8].xyz, frame0.xyz)

        arr0 = traj[6][0]
        for frame0 in traj.iterframe(start=2, step=4, stop=8):
            pass

        for frame0 in traj.iterframe(start=2, step=4, stop=8):
            pass
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        for frame0 in traj.iterframe(start=2, step=2):
            pass
        assert_almost_equal(traj[8].xyz, frame0.xyz)

        count = 0
        for frame0 in traj.iterframe(start=2):
            count += 1
        assert_almost_equal(traj[-1].xyz, frame0.xyz)

        count = 0
        for frame0 in traj.iterframe(start=2, stop=7):
            count += 1
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        for frame0 in traj.iterframe():
            pass
        assert_almost_equal(traj[-1].xyz, frame0.xyz)

        for frame0 in farray.iterframe():
            pass
        assert_almost_equal(traj[-1].xyz, frame0.xyz)

        for frame0 in traj():
            pass
        assert_almost_equal(traj[-1].xyz, frame0.xyz)

        for frame0 in farray():
            pass
        assert_almost_equal(farray[-1].xyz, frame0.xyz)

        count = 0
        for frame0 in traj(start=2, stop=7):
            count += 1
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        count = 0
        for frame0 in farray(start=2, stop=7):
            count += 1
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        count = 0
        for frame0 in farray(2, 7, 1):
            count += 1
        assert_almost_equal(traj[6].xyz, frame0.xyz)

        count = 0
        for frame0 in farray(2, 7, 2):
            count += 1
        assert_almost_equal(traj[6].xyz, frame0.xyz)

    def test_1(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        top2 = traj.top.copy()

        for frame in traj(mask='@CA'):
            pass
        top2.strip("!@CA")
        assert frame.n_atoms == top2.n_atoms

        for frame in traj():
            pass
        assert frame.n_atoms == traj[0].n_atoms

        for frame in traj[:](mask='@CA'):
            pass

        f0 = traj[-1]
        f0.strip(traj.top('!@CA'))
        assert_almost_equal(f0.xyz, frame.xyz)
        assert frame.n_atoms == top2.n_atoms

        for frame in traj[:]():
            pass
        assert frame.n_atoms == traj[0].n_atoms

    def test_IterWithMask(self):
        traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        farray = traj[:]
        for frame in farray(mask='@CA'):
            pass
        f0 = traj[-1]
        f0.strip(traj.top('!@CA'))
        assert_almost_equal(f0.xyz, frame.xyz)


if __name__ == "__main__":
    unittest.main()
