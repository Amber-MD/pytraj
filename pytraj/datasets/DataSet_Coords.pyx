# distutils: language = c++

from pytraj._utils cimport get_positive_idx

cdef class DataSet_Coords(DataSet_1D):
    def __cinit__(self):
        # abstract class, dont' create new object here
        #pass
        # make sure that two pointers pointing to the same address
        # behave like `FrameArray`
        self.baseptr0 = <_DataSet*> self.baseptr_2
        self.baseptr_1 = <_DataSet_1D*> self.baseptr_2
        self._top = Topology()

    def __dealloc__(self):
        # abstract class
        pass

    def __getitem__(self, idxs):
        # TODO : same as Trajin class
        # should combine or inherit or ?
        # return either a Frame instance or FrameArray instance
        
        cdef Frame frame
        cdef FrameArray farray
        cdef int start, stop, step
        cdef int i
        cdef int idx_1

        frame = self.allocate_frame()

        frame.py_free_mem = True

        if self.size == 0:
            raise ValueError("Your FrameArray is empty, how can I index it?")
        if not isinstance(idxs, slice):
            idx_1 = get_positive_idx(idxs, self.size)
            # raise index out of range
            if idxs != 0 and idx_1 == 0:
                # need to check if array has only 1 element. 
                # arr[0] is  arr[-1]
                if idxs != -1:
                    raise ValueError("index is out of range")
            self.baseptr_2.GetFrame(idx_1, frame.thisptr[0])
            return frame
        else:
            # creat a subset array of `FrameArray`
            farray = FrameArray()
            farray.top = self.top
            if idxs.step == None:
                step = 1
            else:
                step = idxs.step
            if idxs.stop == None:
                stop = self.size
            else:
                stop = idxs.stop
            if idxs.start == None:
                start = 0
            else:
                start = idxs.start

            for i in range(start, stop, step):
                # turn `copy` to `False` to have memoryview
                farray.append(self[i], copy=True)
            return farray

    def __setitem__(self, int idx, Frame other):
        idx_1 = get_positive_idx(idx, self.size)
        # raise index out of range
        if idx != 0 and idx_1 == 0:
            # need to check if array has only 1 element. 
            # arr[0] is  arr[-1]
            if idx != -1:
                raise ValueError("index is out of range")
        self.baseptr_2.SetCRD(idx, other.thisptr[0])

    def allocate_frame(self):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.baseptr_2.AllocateFrame()
        return frame

    # do we need this method?
    def _set_topology(self, Topology top):
        self.baseptr_2.SetTopology(top.thisptr[0])

    property top:
        def __get__(self):
            self._top.thisptr[0] = self.baseptr_2.Top()
            return self._top

        def __set__(self, Topology other):
            self.baseptr_2.SetTopology(other.thisptr[0])
    # TODO: add more virtual methods?

    def add_frame(self, Frame frame):
        # TODO : add indices and FrameArray, Trajin_Single ...
        self.baseptr_2.AddFrame(frame.thisptr[0])

    def append(self, frame):
        """alis of addframe"""
        #self.baseptr_2.AddFrame(frame.thisptr[0])
        self.addframe(frame)

    def getframe(self, int idx, Frame frameout):
        self.baseptr_2.GetFrame(idx, frameout.thisptr[0])
