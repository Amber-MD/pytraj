import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        trajlist = []
        for i in range(4):
            traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            trajlist.append(traj)

        mdio.write_traj("./output/test_savedtrajlist.x",
                        traj=trajlist, top=traj.top,
                        overwrite=True)

        # test loading
        traj2 = mdio.load("./output/test_savedtrajlist.x", traj.top)
        N = traj2.size

        ref = traj[0]
        for frame in traj2:
            print (ref.rmsd(frame))

if __name__ == "__main__":
    unittest.main()
