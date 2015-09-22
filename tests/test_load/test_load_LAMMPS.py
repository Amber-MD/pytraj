from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal


class Test(unittest.TestCase):
    def test_0(self):
        # TODO: add assert
        import MDAnalysis as MD
        from MDAnalysis.tests.datafiles import LAMMPSdata
        u = MD.Universe(LAMMPSdata)

        # load to pytraj
        traj = mdio.load_MDAnalysis(u)
        #print(traj)
        #print(traj.top)
        #print(traj.xyz.shape)
        #print(traj.top.atom_names)
        #print(traj.top.residue_names)


if __name__ == "__main__":
    unittest.main()
