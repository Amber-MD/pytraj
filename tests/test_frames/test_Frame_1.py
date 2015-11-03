import unittest
import pytraj as pt
from array import array
import numpy as np
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj import Frame
from pytraj.base import *
from pytraj.math import Vec3
from pytraj import io as mdio

SMALL = 1E-6

N_ATOMS = 10
FRAME = Frame(N_ATOMS)
arr = np.arange(3 * N_ATOMS)
FRAME.set_from_crd(arr)
FRAME_orig = FRAME.copy()


class TestFrame(unittest.TestCase):
    def test_fit(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajnew = mdio.iterload(
            "./data/md1_prod.fit_to_first.Tc5b.x", "./data/Tc5b.top")

        # make sure 0-th frame does not change
        frame0 = traj[0]
        frame0new = trajnew[0]
        assert (frame0.coords == frame0new.coords) == True

        frame1 = traj[1]
        frame1new = trajnew[1]
        assert (frame1.coords == frame1new.coords) == False

        # try do-fitting from Python
        # not right yet
        rmsd, mat, v1, v2 = frame1.rmsd(frame0, get_mvv=True)
        frame1.trans_rot_trans(v1, mat, v2)
        assert frame1.rmsd(frame1new) < 1E-3
        assert frame1new.rmsd(frame1, top=trajnew.top, mask="@CA") < 1E-3

    def test_create_frame(self):
        frame = Frame(list(range(300)))
        assert frame.size == 300
        assert frame.n_atoms == 100
        assert frame.coords == array('d', [x for x in range(300)])
        assert_almost_equal(np.asarray(frame[:, :][0]), frame.coords[:3])

        arr0 = np.arange(300, 0, -1).reshape(100, 3)
        frame[0, :] = array('d', arr0[0, :])
        assert arr0.shape == (100, 3)
        assert frame._buffer2d.shape == (100, 3)
        assert frame[0, 1] == arr0[0, 1]

    def test_buffer1d(self):
        FRAME.buffer1d[0] = 199.
        assert FRAME[0, 0] == FRAME.buffer1d[0] == FRAME.coords[0]

        FRAME[0] = 0.1
        assert FRAME[0, 0] == FRAME.buffer1d[0] == FRAME.coords[0]
        assert FRAME.buffer1d[-1] == FRAME[-1, 2]

        FRAME.buffer1d[:3] = array('d', [1, 2.3, 3.4])
        assert FRAME.coords[:3] == array('d', [1, 2.3, 3.4])

        arr0 = np.asarray(FRAME.buffer1d)
        arr1 = arr0.reshape(10, 3)
        arr1[1] = [100., 200., 300.]
        start = 0
        stop = 10
        strip = 2
        arr = np.asarray(FRAME[start:stop:strip])

    def test_indexing(self):
        # create a Frame instance with N_ATOMS atoms
        N_ATOMS = 10
        frame = Frame(N_ATOMS)
        arr = np.arange(3 * N_ATOMS, dtype=float)
        frame.set_from_crd(arr)

        assert_almost_equal(frame[-2], frame[N_ATOMS - 2])
        assert_almost_equal(frame[-2], arr.reshape(N_ATOMS, 3)[-2])

        frame[-1] = [100., 0, 0]
        assert frame[-1, 0] == 100.

        #frame[-2] = 101.
        #assert frame[-2] == frame[N_ATOMS*3 - 2] == 101.

        frame[0, 0] = 100.
        assert frame[0, 0] == 100.

    def test_other_stuff(self):

        farray = Trajectory(
            "./data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            indices=(1, ))
        frame0 = farray[0]
        atm = AtomMask("@CA")
        farray.top.set_integer_mask(atm)
        frame1 = Frame(frame0, atm)
        frame2 = frame0.get_subframe(mask="@CA", top=farray.top)
        frame3 = frame0.get_subframe("@CA", farray.top)
        assert frame1.coords == frame2.coords == frame3.coords
        frame4 = frame0.get_subframe("!@CA", farray.top)

    def test_rmsd_return_mat_vec_vec(self):
        # TODO : add assert
        farray = Trajectory(
            "./data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            indices=(0, 1))
        frame0 = farray[0]
        rmsd, mat, v1, v2 = frame0.rmsd(farray[1], get_mvv=True)
        assert abs(rmsd - 10.3964) < 1E-3
        arr1 = np.asarray(frame0.buffer1d)[:3]
        frame0.translate(v1)

    def test_long(self):
        N_ATOMS = 10
        # create frame instance with 10 atoms
        frame = Frame(N_ATOMS)
        frameref = Frame(N_ATOMS)

        arr = np.random.rand(N_ATOMS * 3)
        arr_reshape = arr.reshape(N_ATOMS, 3)
        frame.set_from_crd(arr, 30, 0, False)
        assert frame.n_atoms == N_ATOMS
        assert frame.size == N_ATOMS * 3

        assert_almost_equal(np.array(frame.atoms(0)), arr_reshape[0])
        assert_almost_equal(frame[0], arr[:3])

        # frame.info('frame info')
        frame.swap_atoms(1, 8)


        frame.update_atom(1, array('d', [1., 1000., 3000.]))
        frame[3] = 1000000.

        frame.divide(2.)

        count = 0
        for x in frame:
            count += 1
        assert count == N_ATOMS

        for i, x in enumerate(frame):
            if i == 9:
                old_i = frame[i]
                frame[i] = array('d', [1010., 0., 0.])
        assert i == N_ATOMS - 1
        assert_almost_equal(x, old_i)
        assert frame[9][0] == 1010.

        frame.zero_coords()
        arr = np.asarray(frame.coords)
        frame[0] = 1001.10
        assert frame[0, 0] == frame.coords[0]

        arrref = np.random.rand(30)
        frameref.set_from_crd(arr, 30, 0, False)

        frame.update_atoms(
            array('i', [0, 3]), array('d', [0., 0., 0.1, 1.1, 2.3, 3.]))
        assert frame.atoms(0) == array('d', [0., 0., 0.1])
        assert frame.atoms(3) == array('d', [1.1, 2.3, 3.])

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

    def test_tranlate(self):
        farray = Trajectory(
            "./data/md1_prod.Tc5b.x", "./data/Tc5b.top",
            indices=(0, 1))
        f0 = farray[0]
        f1 = f0.copy()
        f2 = f0.copy()

        mylist = [1., 2., 3.]
        vec3 = Vec3(mylist)
        assert isinstance(vec3, Vec3) == True
        assert vec3.tolist() == mylist
        f0.translate(vec3)
        f1.translate(mylist)
        f2.translate(vec3.to_ndarray())
        assert_almost_equal(f0.coords, f1.coords)
        assert_almost_equal(f0.coords, f2.coords)

    def test_velocity_and_force(self):
        traj = pt.load_sample_data('ala3')
        assert not traj[0].has_force(), 'does not have force'
        assert not traj[0].has_velocity(), 'does not have force'
        assert traj[0].force is None, 'force must be None'
        assert traj[0].velocity is None, 'velocity is None'


if __name__ == "__main__":
    unittest.main()
