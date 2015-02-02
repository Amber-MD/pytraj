import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # FrameArray()
        fa = traj[:]
        fname = "./output/test_savemethod.x"
        fname2 = "./output/test_savemethod_2.x"
        fa.save(fname, overwrite=True)
        traj.save(fname2, overwrite=True)

        # load
        fanew = mdio.load(fname, fa.top)
        fanew2 = mdio.load(fname2, fa.top)
        assert fanew.size == fa.size == fanew2.size

        for idx, f0 in enumerate(fa):
            f0new = fanew[idx]
            f0new2 = fanew2[idx]
            assert_almost_equal(f0.coords, f0new.coords)
            assert_almost_equal(f0.coords, f0new2.coords)

    def test_fancy_save(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj[1:8].save("./output/test_fancy_save_frame1_7.x", overwrite=True)

        fanew = mdio.load("./output/test_fancy_save_frame1_7.x", traj.top)

        for idx, f0 in enumerate(traj[1:8]):
            f0new = fanew[idx]
            assert_almost_equal(f0.coords, f0new.coords)

if __name__ == "__main__":
    unittest.main()
