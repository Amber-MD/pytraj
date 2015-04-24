import unittest
from array import array
import numpy as np
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.Frame import Frame
from pytraj.base import *
from pytraj import Vec3
from pytraj import io as mdio
from pytraj.decorators import no_test
from rmsd import rmsd as arr_rmsd

SMALL = 1E-6

N_ATOMS = 10
FRAME = Frame(N_ATOMS)
arr = np.arange(3 * N_ATOMS)
FRAME.set_from_crd(arr)
FRAME_orig = FRAME.copy()

class TestFrame(unittest.TestCase):
    #@no_test
    def test_fit(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        trajnew = mdio.load("./data/md1_prod.fit_to_first.Tc5b.x", "./data/Tc5b.top")

        # make sure 0-th frame does not change
        frame0 = traj[0]
        frame0new = trajnew[0]
        assert (frame0.coords == frame0new.coords) == True
        print(frame0[:10])
        print(frame0new[:10])

        frame1 = traj[1]
        frame1new = trajnew[1]
        print(frame1[:10])
        print(frame1new[:10])
        assert (frame1.coords == frame1new.coords) == False

        # try do-fitting from Python
        # not right yet
        rmsd, mat, v1, v2 = frame1.rmsd(frame0, get_mvv=True)
        print(rmsd)
        frame1.trans_rot_trans(v1, mat, v2)
        #frame1.trans_rot_trans(v1, mat, Vec3(0, 0, 0))
        #assert frame1[:10] == frame1new[:10]
        print(frame1.rmsd_nofit(frame1new))
        print(arr_rmsd(np.asarray(frame1.coords), np.asarray(frame1new.coords)))
        print(frame1.rmsd(frame1new))
        assert frame1.rmsd(frame1new) < 1E-3
        assert frame1new.rmsd(frame1, top=trajnew.top, mask="@CA") < 1E-3

    #@no_test
    def test_strip_atoms(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        frame0 = traj[0]
        frame1 = traj[1]
        frame2 = traj[2]
        frame3 = traj[3]
        frame1 = frame0.strip_atoms("!@CA", traj.top, copy=True)
        frame2 = Frame(frame0)
        print(frame2.n_atoms)
        print(frame1.n_atoms)
        print(frame1)
        assert frame1.n_atoms == 20
        assert frame0.n_atoms == 304
        frame0.strip_atoms("!@CA", traj.top, copy=False)
        assert frame0.n_atoms == 20

        _, mat, v1, v2 = frame2.rmsd(frame3, get_mvv=True)
        print(mat, v1, v2)
        print(frame2[:10])
        print(mdio.writetraj("./output/test_0_before.pdb", traj=frame2, top=traj.top, overwrite=True))
        frame2.trans_rot_trans(v1, mat, v2)
        print(frame2[:10])
        print(mdio.writetraj("./output/test_0_afeter.pdb", traj=frame2, top=traj.top, overwrite=True))
        print(Trajout().help())

    #@no_test
    def test_create_frame(self):
        print("test creating Frame from a list!")
        frame = Frame(list(range(300)))
        assert frame.size == 300
        assert frame.n_atoms == 100
        assert frame.coords == array('d', [x for x in range(300)])
        assert_almost_equal(np.asarray(frame[:, :][0]), frame.coords[:3])

        arr0 = np.arange(300, 0, -1).reshape(100, 3)
        frame[0, :] = array('d', arr0[0, :])
        assert arr0.shape == (100, 3)
        assert frame.buffer2d.shape == (100, 3)
        assert frame[0, 1] == arr0[0, 1]
        ##frame[:], np.arange(300, 0, -1).reshape(100, 3))
        #assert frame.coords == array('d', range(300, 0, -1))

        #frame[:] = np.arange(300, dtype='d')
        #assert frame.coords == array('d', [x for x in range(300)])
        #assert frame[:] == frame.coords
  
    #@no_test
    def test_buffer1d(self):
        print("+++++start test_buffer1d+++++++")
        print(FRAME.coords)
        print(FRAME.buffer1d)

        FRAME.buffer1d[0] = 199.
        print(FRAME.coords[0])
        print(FRAME[0])
        assert FRAME[0, 0] == FRAME.buffer1d[0] == FRAME.coords[0]

        FRAME[0] = 0.1
        assert FRAME[0, 0] == FRAME.buffer1d[0] == FRAME.coords[0]
        assert FRAME.buffer1d[-1] == FRAME[-1, 2]

        FRAME.buffer1d[:3] = array('d', [1, 2.3, 3.4])
        assert FRAME.coords[:3] == array('d', [1, 2.3, 3.4])

        print("test slices")
        print(FRAME[0:10:2]) 

        arr0 = np.asarray(FRAME.buffer1d)
        print(arr0.__len__())
        arr1 = arr0.reshape(10, 3)
        arr1[1] = [100., 200., 300.]
        print(FRAME.atoms)
        #arr0[10] = [100., 200., 300.]
        #assert frame.atoms(10) == arr0[10]

        print("test memoryview for slices")
        # TODO : add
        start = 0
        stop = 10
        strip = 2
        arr = np.asarray(FRAME[start:stop:strip]) 

        # does not work, got "ValueError: ndarray is not contiguous"
        #FRAME[start:stop:strip] = arr_tmp[start:stop:strip]

    #@no_test

    #@no_test
    def test_indexing(self):
        # create a Frame instance with N_ATOMS atoms
        N_ATOMS = 10
        frame = Frame(N_ATOMS)
        arr = np.arange(3 * N_ATOMS, dtype=float)
        frame.set_from_crd(arr)
        print(frame.coords)

        #print "test negative indexing"
        print(frame[-1])
        #assert frame[-1] == 29.
        print(frame[-2])
        assert_almost_equal (frame[-2], frame[N_ATOMS - 2])
        print(np.asarray(frame[-2]))
        print(arr.reshape(N_ATOMS, 3)[-2])
        assert_almost_equal(frame[-2], arr.reshape(N_ATOMS, 3)[-2])

        frame[-1] = [100., 0, 0]
        assert frame[-1, 0] == 100.
        #print frame[-1]

        #frame[-2] = 101.
        #assert frame[-2] == frame[N_ATOMS*3 - 2] == 101.

        print(frame[0:10])
        frame[0, 0] = 100.
        assert frame[0, 0] == 100.

    #@no_test
    def test_other_stuff(self):
        print("print FRAME")
        print(FRAME)

        farray = FrameArray("./data/md1_prod.Tc5b.x", "./data/Tc5b.top", indices=(1,))
        frame0 = farray[0]
        print(frame0)
        atm = AtomMask("@CA")
        farray.top.set_integer_mask(atm)
        print(dir(atm))
        frame1 = Frame(frame0, atm)
        frame2 = frame0.get_subframe(mask="@CA", top=farray.top)
        frame3 = frame0.get_subframe("@CA", farray.top)
        assert frame1.coords == frame2.coords == frame3.coords
        frame4 = frame0.get_subframe("!@CA", farray.top)
        print(frame4)

    #@no_test
    def test_rmsd_return_mat_vec_vec(self):
        # TODO : add assert
        print("test_rmsd_return_mat_vec_vec")
        farray = FrameArray("./data/md1_prod.Tc5b.x", "./data/Tc5b.top", indices=(0,1))
        frame0 = farray[0]
        rmsd, mat, v1, v2 = frame0.rmsd(farray[1], get_mvv=True)
        print(rmsd,  mat, v1, v2)
        assert abs(rmsd - 10.3964) < 1E-3
        arr1 = np.asarray(frame0.buffer1d)[:3]
        print(frame0.coords[:3])
        print(arr1)
        frame0.translate(v1)
        print(v1.tolist())
        print(np.array(v1.tolist(), dtype='d'))
        print(arr1)
        print(frame0[:3])

    #@no_test
    def test_long(self):
        N_ATOMS = 10
        # create frame instance with 10 atoms
        frame = Frame(N_ATOMS)
        frameref = Frame(N_ATOMS)
        
        arr = np.random.rand(N_ATOMS * 3)
        arr_reshape = arr.reshape(N_ATOMS, 3)
        frame.set_from_crd(arr, 30, 0, False)
        assert frame.n_atoms == N_ATOMS
        assert frame.size == N_ATOMS*3
        
        print(frame.atoms(0))
        print(arr_reshape[0])
        assert_almost_equal(np.array(frame.atoms(0)), arr_reshape[0])
        assert_almost_equal(frame[0], arr[:3])

        # frame.info('frame info')
        frame.swap_atoms(1, 8)
        
        print("after swapping 1 - 8")

        print("update coords_copy for atom 1")
        frame.update_atom(1, array('d', [1., 1000., 3000.]))
        print(frame.atoms(1))
        print(frame[3])
        print("assign frame[3] to 1000000.")
        frame[3] = 1000000.
        print(frame.atoms(1))
        
        print("deviding Frame")
        frame.divide(2.)

        print("test iteration")
        count = 0
        for x in frame:
            count += 1
        print("count = %s " % count)
        assert count == N_ATOMS

        print("test enumeration")
        for i, x in enumerate(frame):
            if i == 9:
                old_i = frame[i]
                frame[i] = array('d', [1010., 0., 0.])
        print(i)
        assert i == N_ATOMS - 1
        print(x)
        assert_almost_equal (x, old_i)
        assert frame[9][0] == 1010.

        print("set zero_coords_copy")
        frame.zero_coords()
        arr = np.asarray(frame.coords)
        frame[0] = 1001.10
        assert frame[0, 0] == frame.coords[0]
        print(frame[0])
        
        arrref = np.random.rand(30)
        frameref.set_from_crd(arr, 30, 0, False)

        frame.update_atoms(array('i', [0, 3]), array('d', [0., 0., 0.1, 1.1, 2.3, 3.]))
        print(type(frame.atoms(0)))
        print(frame.atoms(0))
        assert frame.atoms(0) == array('d', [0., 0., 0.1])
        assert frame.atoms(3) ==  array('d', [1.1, 2.3, 3.])

    def test_iter(self):
        print("test iteration")
        alist = []
        frame = FRAME_orig.copy()

        for x in frame:
            alist += [int(a) for a in x]
        print(alist)
        assert alist == list(range(3 * N_ATOMS))

        print("test enumerate")
        alist = []
        for idx, x in enumerate(frame):
            alist += [int(a) for a in x]
        assert alist == list(range(3 * N_ATOMS))
        print(alist)
        print("====================end test_iter")

if __name__ == "__main__":
    unittest.main()
