# distutils: language = c++
from pytraj.decorators import for_testing


cdef class DataSet_Coords_TRJ(DataSet_Coords):
    def __cinit__(self, *args):
        # TODO : do we really need to epoxe _DataSet and _DataSet_1D?
        # seriouly I need to use 4 pointers for class inheritance.
        # use pointer casting instead? (look ugly?)
        self.baseptr0 = <_DataSet*> new _DataSet_Coords_TRJ()
        # recast
        self.baseptr_1 = <_DataSet_1D*> self.baseptr0
        self.baseptr_2 = <_DataSet_Coords*> self.baseptr0
        self.thisptr = <_DataSet_Coords_TRJ*> self.baseptr0

    def __dealloc__(self):
        del self.thisptr
    
    def __iter__(DataSet_Coords_TRJ self):
        """iterately getting Frame instance
        TODO : get memoery view or copy?
        """
        cdef int i 
        cdef Frame frame
        frame = self.allocate_frame()

        for i in range(self.size):
            self.thisptr._GetFrame(i, frame.thisptr[0])
            yield frame

    def _recast(self):
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_1 = <_DataSet_1D*> self.thisptr
        self.baseptr_2 = <_DataSet_Coords*> self.thisptr

    @classmethod
    def alloc(cls):
        """return base class: DataSet"""
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_Coords_TRJ.Alloc()
        return dset

    def load(self, filename, Topology top=Topology(), ArgList arglist=ArgList()):
        if top.is_empty():
            if not self.top.is_empty():
                top = self.top
            else:
                raise ValueError("need to have non-empty topology file")

        filename = filename.encode())
        if self.top.is_empty() and not top.is_empty():
            print "assigning new non-empty Topology"
            self.top = top
        return self.thisptr.AddSingleTrajin(filename, arglist.thisptr[0], top.thisptr)

    def add_trajin(self, Trajin trajin):
        """add memoryview for input trajin"""
        self.thisptr.AddInputTraj(trajin.baseptr_1)

    @property
    def size(self):
        return self.thisptr.Size()

    def get_frame(self, int idx, Frame frame_in, *args):
        cdef AtomMask atm_in
        if self.top.n_atoms != frame_in.n_atoms:
            raise ValueError("n_atoms should be matched between Frame and Topology")
        if not args:
            self.thisptr._GetFrame(idx, frame_in.thisptr[0])
        else:
            atm_in = args[0]
            self.thisptr._GetFrame(idx, frame_in.thisptr[0], atm_in.thisptr[0])
