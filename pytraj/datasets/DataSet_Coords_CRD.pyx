# distutils: language = c++

from pytraj.Trajin_Single cimport Trajin_Single
from pytraj.Topology cimport Topology
from pytraj.externals.six import string_types


cdef class DataSet_Coords_CRD (DataSet_Coords):
    def __cinit__(self):
        self.thisptr = new _DataSet_Coords_CRD()
        self.baseptr0 = <_DataSet*> self.thisptr
        self.baseptr_1 = <_DataSet_Coords*> self.thisptr

        # let python frees memory
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    @classmethod
    def alloc(self):
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_Coords_CRD.Alloc()
        return dset

    def load(self, filename_or_traj, top=Topology(), copy_top=False):
        cdef Topology tmp_top
        cdef Trajin_Single trajin_single
        cdef Frame frame

        if isinstance(top, string_types):
            top = Topology(top)

        if top.is_empty():
            if not self.top.is_empty():
                tmp_top = self.top
            else:
                raise ValueError("need to have non-empty topology file")
        else:
            tmp_top = top
            # update self.top too
            if copy_top == True:
                self.top = top.copy()
            else:
                self.top = top

        if isinstance(filename_or_traj, string_types):
            trajin_single = Trajin_Single(filename_or_traj, tmp_top)
            for frame in trajin_single:
                self.append(frame)
        else:
            # assume that we can iterate over filename_or_traj to get Frame object
            for frame in filename_or_traj:
                self.append(frame)
