import unittest
import numpy as np
import pytraj as pt
from pytraj import Frame
from pytraj.core import Box
from pytraj.testing import aa_eq
from pytraj.compat import zip


class TestBox(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame0.box.tolist()
        frame0.box = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 6.])
        assert frame0.has_box(), 'must has box'
        frame0.set_nobox()
        assert not frame0.has_box(), 'must not has box when setting no box'

    def test_1(self):
        box = Box()
        box.type = 'truncoct'
        box.set_nobox()

        dummy = 100.
        box.data[0] = dummy
        assert box.data[0] == dummy
        assert box.values[0] == dummy

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
        aa_eq(f1.box.values, box.values, decimal=7)

    def test_real_box(self):
        traj = pt.load("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        trajiter = pt.iterload("./data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        saved_box = Box([3.94559740E+01, 4.68215170E+01, 4.04695410E+01, 90.,
                         90., 90.])
        aa_eq(traj.top.box.values, saved_box.values)
        for frame in traj:
            assert frame.box.type == 'ortho'
            aa_eq(frame.box.values,
                  [
                      35.2627796623, 41.8455476799, 36.168629529, 90.0, 90.0,
                      90.0
                  ],
                  decimal=1)

        arr0 = traj.unitcells
        arr1 = trajiter.unitcells

        for b0, b1, frame in zip(arr0, arr1, trajiter):
            box = frame.box
            # FIXME:
            # python3)
            b2 = box.values
            aa_eq(b0, b1)
            aa_eq(b0, b2)

        # test assign Box with list/tuple
        b = Box(saved_box.values)
        b2 = Box((t for t in saved_box.values))
        aa_eq(b.values, saved_box.values, decimal=7)
        aa_eq(b2.values, saved_box.values, decimal=7)
        # assign frame.box with list/tuple
        frame.box = b.values
        b3 = frame.box
        aa_eq(b3.values, saved_box.values, decimal=7)

    def test_nobox(self):
        from pytraj import Trajectory
        traj = Trajectory()
        traj._allocate(10, 10)

    def test_from_matrix_3x3(self):
        from pytraj.math import Matrix_3x3
        mat = Matrix_3x3()
        box = Box(mat)

    def test_box_constructor_with_type(self):
        # TODO: assertion
        box = Box('rhombic')
        box = Box('ortho')


if __name__ == "__main__":
    unittest.main()
