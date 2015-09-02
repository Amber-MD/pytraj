# distutils: language = c++

from .._shared_methods import _frame_iter
from .._shared_methods import _xyz, _tolist
from .._shared_methods import my_str_method
from ..utils import _import_numpy, _import

from .._cyutils import get_positive_idx

cdef class DataSet_Coords(DataSet):
    def __cinit__(self):
        # abstract class, dont' create new object here
        #pass
        # make sure that two pointers pointing to the same address
        self.baseptr0 = <_DataSet*> self.baseptr_1
        self._top = Topology()

    def __dealloc__(self):
        # abstract class
        pass

    @property
    def n_frames(self):
        return self.size

    @property
    def n_atoms(self):
        """used for frame_iter"""
        return self.top.n_atoms

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args, **kwd):
        return self.frame_iter(*args, **kwd)

    def __iter__(self):
        """iterately getting Frame instance
        TODO : get memoryview or copy?
        """
        cdef int i 
        cdef Frame frame
        frame = self.allocate_frame()

        for i in range(self.size):
            self.baseptr_1.GetFrame(i, frame.thisptr[0])
            yield frame

    def __getitem__(self, idx):
        cdef Frame frame
        frame = self.allocate_frame()
        frame.py_free_mem = True

        if self.size == 0:
            raise ValueError("Your Trajectory is empty, how can I index it?")
        self.baseptr_1.GetFrame(idx, frame.thisptr[0])
        self.tmpfarray = frame
        return self.tmpfarray

    def __setitem__(self, int idx, Frame other):
        idx_1 = get_positive_idx(idx, self.size)
        # raise index out of range
        if idx != 0 and idx_1 == 0:
            # need to check if array has only 1 element. 
            # arr[0] is  arr[-1]
            if idx != -1:
                raise ValueError("index is out of range")
        self.baseptr_1.SetCRD(idx, other.thisptr[0])

    def frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
        return _frame_iter(self, start, stop, stride, mask)

    def allocate_frame(self):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.baseptr_1.AllocateFrame()
        return frame

    property top:
        def __get__(self):
            self._top.thisptr[0] = self.baseptr_1.Top()
            return self._top

        def __set__(self, Topology other):
            self.baseptr_1.SetTopology(other.thisptr[0])

    def add_frame(self, Frame frame):
        self.baseptr_1.AddFrame(frame.thisptr[0])

    def append(self, frame):
        """alis of addframe"""
        self.add_frame(frame)

    def get_frame(self, int idx, Frame frameout):
        self.baseptr_1.GetFrame(idx, frameout.thisptr[0])

    @property
    def xyz(self):
        """return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object
        """
        cdef Frame frame
        cdef int i
        _, np = _import_numpy()
        n_frames = self.n_frames 
        n_atoms = self.top.n_atoms
        arr = np.empty((n_frames, n_atoms, 3))

        for i in range(n_frames):
            arr[i] = self[i].xyz
        return arr

    def tolist(self):
        """return flatten list for traj-like object"""
        cdef Frame frame
        return [frame.tolist() for frame in self]

    def to_dataframe(self):
        raise NotImplementedError()
