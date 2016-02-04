import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):

    def test_0(self):
        trajlist = []
        N = 4
        for i in range(N):
            traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
            trajlist.append(traj)

        traj0 = trajlist[0]

        mdio.write_traj("./output/test_savedtrajlist.x",
                        traj=trajlist,
                        top=traj.top,
                        overwrite=True)

        # test loading
        traj2 = mdio.iterload("./output/test_savedtrajlist.x", traj.top)
        Nsize = int(traj2.n_frames / 4)
        traj0_new = traj2[:Nsize]

        for frame0, frame0_new in zip(traj0, traj0_new):
            assert (frame0.rmsd(frame0_new) < 1E-3)


if __name__ == "__main__":
    unittest.main()
