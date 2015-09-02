# distutils: language = c++

from ..trajs.TrajectoryCpptraj import TrajectoryCpptraj
from ..Topology cimport Topology
from ..externals.six import string_types

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

    def load(self, filename_or_traj, top=Topology(), copy_top=False, copy=True):
        cdef Topology tmp_top
        cdef Frame frame

        if isinstance(top, string_types):
            self.top = top = Topology(top)

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
            trajin_single = TrajectoryCpptraj()
            trajin_single.load(filename_or_traj, tmp_top)
            for frame in trajin_single:
                self.append(frame.copy()) # always copy
        else:
            # assume that we can iterate over filename_or_traj to get Frame object
            for frame in filename_or_traj:
                if copy:
                    self.append(frame.copy())
                else:
                    self.append(frame)
