from __future__ import print_function
import unittest
from pytraj.decorators import no_test
from pytraj.base import *
from pytraj import adict, cast_dataset
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.six_2 import izip

class Test(unittest.TestCase):
    @no_test
    def test_0(self):
        print ("coords")
        dslist = DataSetList()
        coords = cast_dataset(dslist.add_set("coords", "my_coords", ""), dtype='coords')
        trajin_and_top = ("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.load(*trajin_and_top)
        coords.load(*trajin_and_top)
        assert coords.size == traj.size

        for f0, f1 in izip(coords(2, 8, 1), traj(2, 8, 1)):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

        for f0, f1 in izip(coords(), traj()):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

    def test_1(self):
        print ("coords_trj")
        dslist = DataSetList()
        dslist.add_set("traj", "my_coords", "")
        trajin_and_top = ("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj = mdio.load(*trajin_and_top)
        coords = dslist[0]
        coords.load(*trajin_and_top)
        print (coords.top)
        print (coords)
        assert coords.size == traj.size

        for i in range(coords.size):
            print (coords[i])

        for f0, f1 in izip(coords(2, 8, 1), traj(2, 8, 1)):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

        for f0, f1 in izip(coords(), traj()):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

        # FIXME: frame_iter for TRJ and CRD results infinite loop
        # ack
        for f0, f1 in izip(coords.frame_iter(2, 8, 1, "!@CA"), traj(2, 8, 1, '!@CA')):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

        for f0, f1 in izip(coords(mask='!@CA'), traj(mask='!@CA')):
            print (f0[0, :], f1[0, :])
            assert_almost_equal(f0.coords, f1.coords)

if __name__ == "__main__":
    unittest.main()
