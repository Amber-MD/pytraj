import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.core import Box
from array import array as pyarray
from pytraj.decorators import test_if_having
from pytraj.testing import eq, aa_eq
from pytraj.compat import zip


class TestBox(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0.box_crd()
        frame0.boxview[:] = pyarray('d', [0.0, 1.0, 2.0, 3.0, 4.0, 6.])
        frame0.set_nobox()

    def test_1(self):
        box = Box()
        box.set_trunc_oct()
        box.set_nobox()

        dummy = 100.
        box.data[0] = dummy
        assert box.data[0] == dummy
        assert box.tolist()[0] == dummy

    @test_if_having("numpy")
    def test_2(self):
        import numpy as np
        box = Box()
        arr0 = np.arange(6).astype(np.float64)
        box.data[:] = arr0

        for idx, x in enumerate(arr0):
            assert box.data[idx] == x

        # set Box for Frame
        f1 = Frame()
        f1.box = box
        eq(f1.box.tolist(), box.tolist())

        f2 = Frame()
        f2.box = box.to_ndarray()
        eq(f2.box.tolist(), box.tolist())

    def test_real_box(self):
        traj = mdio.load("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        trajiter = mdio.iterload("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        saved_box = Box(
            [3.94559740E+01, 4.68215170E+01, 4.04695410E+01, 90., 90., 90.])
        #print(traj.top.box)
        #print(trajiter.top.box)
        aa_eq(traj.top.box.tolist(), saved_box.tolist())
        for frame in traj:
            assert frame.box.type == 'ortho'
            aa_eq(frame.box.tolist(), [
                35.2627796623, 41.8455476799, 36.168629529, 90.0, 90.0, 90.0
            ],
                  decimal=1)

        arr0 = traj.unitcells
        arr1 = trajiter.unitcells

        for b0, b1, frame in zip(arr0, arr1, trajiter):
            box = frame.box
            # FIXME:
            # b2 = frame.box.to_ndarray() # got wrong box in python2 (ok with
            # python3)
            b2 = box.to_ndarray()
            aa_eq(b0, b1)
            aa_eq(b0, b2)

        # test assign Box with list/tuple
        b = Box(saved_box.tolist())
        b2 = Box((t for t in saved_box.tolist()))
        assert (b.tolist() == saved_box.tolist()) == True
        assert (b2.tolist() == saved_box.tolist()) == True
        # assign frame.box with list/tuple
        frame.box = b.tolist()
        b3 = frame.box
        assert (b3.tolist() == saved_box.tolist()) == True

    def test_nobox(self):
        from pytraj import Trajectory
        traj = Trajectory()
        traj._allocate(10, 10)
        #print(traj.unitcells)

    def test_assign_box_type(self):
        #print("test_assign_box_type")
        box = Box()
        assert box.type == 'nobox'
        box.type = 'ortho'
        assert box.type == 'ortho'
        assert box.alpha > 0.
        box.type = 'truncoct'
        assert box.type == 'truncoct'

        box.type = 'rhombic'
        assert box.type == 'rhombic'

        #print("rhombic box?")
        #print(box)

        # assert raise if not correctly set type
        def wrong_word(box=box):
            box.type = 'test'

        self.assertRaises(ValueError, lambda: wrong_word())

        # test update boxtype
        box = Box()
        box.values[3:] = [90., 90., 90.]
        assert box.type == 'nobox'
        box.update_box_type()
        assert box.type == 'ortho'

    def test_from_matrix_3x3(self):
        from pytraj.math import Matrix_3x3
        mat = Matrix_3x3()
        box = Box(mat)

    def test_set_box_from_array(self):
        box = Box()
        box.set_box_from_array([30., 30., 30., 90., 90., 90.])
        assert box.type == 'ortho'

        box.set_box_from_array([30., 30., 30., 60., 90., 60.])
        assert box.type == 'rhombic'

    def test_box_constructor_with_type(self):
        box = Box('rhombic')
        assert box.type == 'rhombic'
        box = Box('ortho')
        assert box.type == 'ortho'


if __name__ == "__main__":
    unittest.main()
