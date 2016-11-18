import unittest

from pytraj import io as pt
from pytraj.testing import tempfolder

from utils import fn


class TestTrajlist(unittest.TestCase):

    def test_trajlist(self):
        trajlist = []
        N = 4
        for i in range(N):
            traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
            trajlist.append(traj)

        traj0 = trajlist[0]

        with tempfolder():
            pt.write_traj("test_savedtrajlist.x",
                            traj=trajlist,
                            top=traj.top,
                            overwrite=True)

            # test loading
            traj2 = pt.iterload("test_savedtrajlist.x", traj.top)
            Nsize = int(traj2.n_frames / 4)
            traj0_new = traj2[:Nsize]

            for frame0, frame0_new in zip(traj0, traj0_new):
                assert (frame0.rmsd(frame0_new) < 1E-3)


if __name__ == "__main__":
    unittest.main()
