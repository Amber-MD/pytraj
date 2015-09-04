import unittest
from array import array
from pytraj.base import *
from pytraj import io as mdio
import numpy as np
from pytraj.utils import assert_almost_equal


class TestCHARMM(unittest.TestCase):
    def test_0(self):
        top = Topology("./data/ala3.psf")
        top.summary()
        reslit = top.residuelist
        atm = AtomMask("@CA")
        top.set_integer_mask(atm)

        atm.invert_mask()
        frame = Frame(atm.n_atoms)
        frame[:10] = np.asarray(array('d', list(range(30)))).reshape(10, 3)
        frame.zero_coords()

    def test_1(self):
        traj = mdio.iterload("./data/ala3.dcd", "./data/ala3.psf")
        traj.save("./output/_save_charmm_to_amber.x", overwrite=True)
        # test loading
        trajamber = mdio.iterload("./output/_save_charmm_to_amber.x",
                                  "./data/ala3.psf")
        for i in range(traj.n_frames):
            assert_almost_equal(trajamber[i].coords, traj[i].coords)


if __name__ == "__main__":
    unittest.main()
