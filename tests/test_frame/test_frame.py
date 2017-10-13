import unittest
import pytraj as pt
from utils import fn
from array import array
import numpy as np
from pytraj import Frame

from pytraj.testing import aa_eq
from pytraj import *

N_ATOMS = 10
FRAME = Frame(N_ATOMS)
arr = np.arange(3 * N_ATOMS)
FRAME.xyz[:] = arr.reshape(N_ATOMS, 3)
FRAME_orig = FRAME.copy()


class TestFrame(unittest.TestCase):
    def test_fit(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        trajnew = pt.iterload(
            fn('md1_prod.fit_to_first.Tc5b.x'), fn('Tc5b.top'))

        # make sure 0-th frame does not change
        frame0 = traj[0]
        trajnew[0]

        frame1 = traj[1]
        frame1new = trajnew[1]

        # try do-fitting from Python
        # not right yet
        rmsd, mat, v1, v2 = frame1.rmsd(frame0, get_mvv=True)
        frame1._trans_rot_trans(v1, mat, v2)
        assert frame1.rmsd(frame1new) < 1E-3
        assert frame1new.rmsd(frame1, top=trajnew.top, mask="@CA") < 1E-3

    def test_buffer1d(self):
        FRAME._buffer1d[0] = 199.
        assert FRAME[0, 0] == 199.

        FRAME[0] = 0.1
        assert FRAME[0, 0] == 0.1
        assert FRAME._buffer1d[-1] == FRAME[-1, 2]

        FRAME._buffer1d[:3] = np.array([1, 2.3, 3.4], dtype='f8')
        aa_eq(FRAME.xyz[0], array('d', [1, 2.3, 3.4]))

        arr0 = np.asarray(FRAME._buffer1d)
        arr1 = arr0.reshape(10, 3)
        arr1[1] = [100., 200., 300.]
        start = 0
        stop = 10
        strip = 2
        arr = np.asarray(FRAME[start:stop:strip])

    def test_rmsd_return_mat_vec_vec(self):
        # TODO : add assert
        farray = Trajectory(fn('Tc5b.x'), fn('Tc5b.top'))
        frame0 = farray[0]
        rmsd, mat, v1, v2 = frame0.rmsd(farray[1], get_mvv=True)
        assert abs(rmsd - 10.3964) < 1E-3
        arr1 = np.asarray(frame0._buffer1d)[:3]

    def test_iter(self):
        alist = []
        frame = FRAME_orig.copy()

        for x in frame:
            alist += [int(a) for a in x]
        assert alist == list(range(3 * N_ATOMS))

        alist = []
        for idx, x in enumerate(frame):
            alist += [int(a) for a in x]
        assert alist == list(range(3 * N_ATOMS))

    def test_velocity_and_force(self):
        traj = pt.load_sample_data('ala3')
        assert not traj[0].has_force(), 'does not have force'
        assert not traj[0].has_velocity(), 'does not have force'
        assert traj[0].force is None, 'force must be None'
        assert traj[0].velocity is None, 'velocity is None'

    def test_velocity_and_force_allocation(self):
        frame = pt.Frame()
        top = pt.tools.make_fake_topology(100)

        assert frame.n_atoms == 0
        assert not frame.has_force()
        assert not frame.has_velocity()

        crdinfo = dict(has_force=True, has_velocity=True)

        frame._allocate_force_and_velocity(top, crdinfo)
        assert frame.n_atoms == 100

        assert frame.has_force()
        assert frame.has_velocity()

    def test_velocity(self):
        frame = pt.Frame(304)
        velocity = np.zeros(304 * 3) + 1.
        assert frame.velocity is None
        frame.velocity = velocity
        assert frame.velocity.shape == (304, 3)
        frame.velocity = None
        assert frame.velocity is None

    def test_force(self):
        frame = pt.Frame(304)
        force = np.zeros(304 * 3) + 1.
        assert frame.force is None
        frame.force = force
        assert frame.force.shape == (304, 3)
        frame.force = None
        assert frame.force is None


if __name__ == "__main__":
    unittest.main()
