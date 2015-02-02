import unittest 
import numpy as np
from time import time
from pytraj.base import *
from pytraj.FrameArray2 import FrameArray2 as TrajectoryReadOnly
from pytraj.decorators import no_test

TRAJ = TrajectoryReadOnly()
TRAJ.top = Topology("data/Tc5b.top")
TRAJ.load("./data/md1_prod.Tc5b.x")

TRAJ2 = FrameArray2()
TRAJ2.top = TRAJ.top
TRAJ2.load("./data/md1_prod.Tc5b.x")

class TestFrameArray2(unittest.TestCase):
    #@no_test
    def test_size_0(self):
        TRAJ_ = TrajectoryReadOnly()
        assert  TRAJ_.size == 0

    #@no_test
    def test_addframe(self):
        # TODO : does not add any Frame
        # TRAJ2 is still 0
        TRAJ2 = TrajectoryReadOnly()
        TRAJ2.top = TRAJ.top.copy()
        TRAJ2.append(TRAJ[0])
        print(TRAJ2.size)

    #@no_test
    def test_slice(self):
        farray0 = TRAJ[:10]
        farray1 = TRAJ[:10]

        farray0[0][0, 0] = 10.
        print(farray0[0, 0, 0])
        print(farray1[0, 0, 1])
        # make sure that farray0 and farray1 are NOT views of TRAJ[:10]
        assert farray0[0, 0, 0] == 10.
        assert farray0[0, 0, 0] != farray1[0, 0, 0]

    #@no_test
    def test_indices_0(self):
        # FIXME : incorrect memoryview when using TRAJ2[0].buffer1d
        # why?  TRAJ2[0] will return a temp Frame instance, if we don't hold this instance, 
        # its memory will be freed --> zero coords since we no longer have buffer1d?
        # everytime we use TRAJ2[0], Python actually create a new Frame instance
        # Solution: don't let Python free memory if using this style?
        TRAJ2[0].py_free_mem = False
        buffer1d = TRAJ2[0].buffer1d
        print("start testing `buffer1d`+++++++++")
        print(buffer1d)
        print(buffer1d)
        print("end testing `buffer1d`+++++++++")
        print("calling TRAJ[0] several time, what do we expect?")
        print(TRAJ[0])
        print(TRAJ[0])
        frame = TRAJ[0]
        print(TRAJ[0])
        print(frame)
        print("calling TRAJ[0] several time, is what we are expecting?")
        print(np.asarray(frame.buffer1d)[0:10])
        print(np.asarray(TRAJ[0].buffer1d)[0:10])

        print("test buffer1d of Frame")
        frame.py_free_mem = False
        frame.buffer1d
        print("end test buffer1d of Frame")

    #@no_test
    def test_indices_1(self):
        """test repeating indexing
        """
        frame0_0 = TRAJ2[0]
        frame1 = TRAJ2[0]
        frame2 = TRAJ2[0]
        print(frame0_0, frame1, frame2)
        print(frame0_0.coords[:10])

        assert frame0_0.coords == frame1.coords == frame2.coords

        print("updating coords for frame0_0, make sure does not change frame1, frame2")
        frame0_0[0] = 10000.
        assert frame0_0.coords != frame1.coords == frame2.coords
        assert frame1.coords == frame2.coords
        assert frame0_0.coords != TRAJ2[0].coords
        assert frame1.coords == TRAJ2[0].coords == frame2.coords

        frame0_0_arr = np.asarray(frame0_0.buffer1d)
        print(frame0_0_arr[:10])
        print(frame0_0.coords[:10])
        #assert frame0_0_arr[:10] == np.array(frame0[:10])

    #@no_test
    def test_indices(self):
        """play with frame and traj[0] (frame = traj[0])"""
        traj = TRAJ
        print(traj.size)
        frame = traj[0]
        print(frame)
        print(traj[0])
        print(frame.coords[:10])
        print(traj[0].coords[:10])

        print("make sure frame is not an alias of traj[0]")
        assert frame != traj[0]

        # other test
        print(frame.size, traj[0].size)
        assert frame.size == traj[0].size
        assert frame.n_atoms == traj[0].n_atoms 

        print("make sure frame.buffer1d is equal (but not an alias) of traj[0]")
        frame_arr = np.asarray(frame.buffer1d)

        # can not make buffer1d view for traj0_arr since `traj[0]` will create 
        # a temp Frame instance, when this instance is deallocated, the buffer1d is gone
        # traj0_arr = np.asarray((traj[0]).buffer1d)
        print(frame_arr[:10])
        # print traj0_arr[:10]
        #np.testing.assert_almost_equal(frame_arr, traj0_arr, decimal=3)

    #@no_test
    def test_strip_atoms(self):
        print("test_strip_atoms")
        farray = FrameArray()
        #farray.get_frames(TRAJ, update_top=True)
        farray.get_frames(TRAJ, indices=(1, 3, 5, 7,), update_top=True)
        t0 = time()
        print(farray.top.n_atoms)
        farray.strip_atoms(mask=":2-10@CA", update_top=True)
        print(time() - t0)
        print(farray.size)

    #@no_test
    def test_general(self):
        traj = FrameArray2()
        traj2 = FrameArray2()
        print(dir(traj))
        print(traj.top.is_empty())
        traj.top = Topology("data/Tc5b.top")
        traj2.top = Topology("data/Tc5b.top")
        
        traj.load("./data/md1_prod.Tc5b.x")
        
        print(traj.data_format)
        print(traj.data_type)
        print(traj.is_torsion_array())
        print(traj.is_empty())
        print(traj.ndim)
        print(traj.idx)
        
        frame = traj.allocate_frame()
        print(frame.n_atoms)
        print(frame.size)
        print(frame.is_empty())
        print(frame)
        
        #traj.add_frame(frame)
        print(frame[0])
        traj.get_frame(0, frame)
        print(frame[0])
        
        frame0 = traj[0]
        assert frame0.coords == frame.coords
        
        # make sure that frame and frame0 are copies. Changing frame0 does not affact frame
        frame0[10] = 10000.
        assert frame0.coords != frame.coords
        #traj.get_frame(1, frame)
        #print frame[0]

if __name__ == "__main__":
    unittest.main()
