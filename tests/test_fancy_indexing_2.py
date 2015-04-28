import unittest
from pytraj import io as mdio
from pytraj import FrameArray
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        FA = traj[:]

        print (traj['@CA'])
        frame0 = traj[0]
        print (hasattr(frame0, 'shape'))
        aa_eq(frame0[traj.top("@CA")].flatten(), 
              traj['@CA'].xyz.flatten())

        # slicing with list or array
        indices = [1, 2]
        fa = traj[indices]
        fa2 = FA[indices]
        self.assertIsInstance(fa, FrameArray)
        # from TrajReadOnly
        aa_eq(fa[0].coords, traj[1].coords)
        aa_eq(fa[1].coords, traj[2].coords)
        # from FrameArray
        aa_eq(fa2[1].coords, traj[2].coords)
        aa_eq(fa2[0].coords, traj[1].coords)

if __name__ == "__main__":
    unittest.main()
