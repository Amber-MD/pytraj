import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils import assert_almost_equal as aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = mdio.loadpdb_rcsb("2KOC")
        # coords of 1st atom of 1st frame
        # http://www.rcsb.org/pdb/files/2KOC.pdb
        assert traj.n_frames == 20
        aa_eq(traj[0, 0], [-8.886, -5.163, 9.647])


if __name__ == "__main__":
    unittest.main()
