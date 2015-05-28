import unittest
from array import array
from pytraj.base import *
from pytraj import io as mdio
from numpy.testing import assert_almost_equal
import numpy as np

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[9]
        frame = Frame()
        frame.set_frame_from_mask("@CA", traj.top)
        assert frame.size == 60
        assert frame.n_atoms == 20
        print(frame.n_atoms)
        print(frame.size)
        print(frame[:].shape)
        print(frame[0])
        print(dir(frame))
        frame.zero_coords()
        print(frame[0])
        arr0 = np.asarray([[100, 200, 3],], dtype=np.float64)
        frame.append_xyz(arr0)
        print ("frame from test_0")
        assert_almost_equal(frame[-1], arr0[0], decimal=4)
        #print frame.t_address()

    def test_indexing_AtomMask(self):
        print("test_indexing_AtomMask")
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        atm = AtomMask("@CA")
        traj.top.set_integer_mask(atm)
        frame0 = traj[9]
        arr0_0 = frame0[atm]

        # test AtomSelect
        # turn off since we don't use this class anymore
        #asl = AtomSelect(top=traj.top)
        #asl.selected_frame = traj[9]
        #arr0_1 = asl.select("@CA")

        npassert = np.testing.assert_almost_equal
        #npassert(arr0_0, arr0_1, decimal=5)

        # test dict
        arr0_2 = frame0[dict(top=traj.top, mask='@CA')]
        arr0_3 = frame0[{'top':traj.top, 'mask':'@CA'}]
        npassert(arr0_0, arr0_2, decimal=5)
        npassert(arr0_0, arr0_3, decimal=5)

        # 
        print(frame0[AtomMask(303)])
        print(dir(AtomMask))

if __name__ == "__main__":
    unittest.main()
