import unittest; import pytraj as pt
from time import time
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test
from pytraj import adict


class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = Trajectory()
        farray.top = traj.top

        for i in range(10):
            for frame in traj:
                farray.append(frame)
        print(farray)

        farray.join(farray.copy())

    def test_1(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        farray = Trajectory(traj, traj.top)
        t0 = time()

        farray2 = Trajectory()
        farray2.top = traj.top.copy()

        t0 = time()
        for frame in farray:
            farray2.append(frame)

        t0 = time()
        for frame in farray:
            farray2.append(frame)

        farray3 = farray2[:1000]
        farray3.join((farray3[:], farray3[:], farray3[:]))
        t0 = time()

        dslist = DataSetList()
        for f0 in farray3:
            assert f0.n_atoms == farray3.top.n_atoms
        # adict['distance'](":2@CA :10@CA", farray3,
        #                  farray3.top, dslist=dslist)

        d0 = adict['distance'](":2@CA :10@CA", farray3, farray3.top,
                               quick_get=True)


if __name__ == "__main__":
    unittest.main()
