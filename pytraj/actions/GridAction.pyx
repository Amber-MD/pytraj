# distutils: language = c++
from pytraj.cpptraj_dict import GridModeDict, get_key


cdef class GridAction:
    def __cinit__(self):
        self.thisptr = new _GridAction()

    def __dealloc__(self):
        del self.thisptr

    def init(self, char* callingRoutine, ArgList arglist, DataSetList dslist):
        cdef DataSet_GridFlt dset_gf = DataSet_GridFlt()
        # since we use memory view, we let cpptraj free memory
        dset_gf.py_free_mem = False
        dset_gf.thisptr = self.thisptr.GridInit(callingRoutine, arglist.thisptr[0], dslist.thisptr[0])
        return dset_gf

    def info(self, DataSet_GridFlt dset_gf):
        self.thisptr.GridInfo(dset_gf.thisptr[0])

    def process(self, Topology top):
        return self.thisptr.GridSetup(top.thisptr[0])

    def frame(self, Frame frame, AtomMask atm, DataSet_GridFlt dset_gf):
        self.thisptr.GridFrame(frame.thisptr[0], atm.thisptr[0], dset_gf.thisptr[0])

    def mode(self):
        return get_key(self.thisptr.GridMode(), GridModeDict)

    def center_mask(self):
        cdef AtomMask atm = AtomMask()
        atm.thisptr[0] = self.thisptr.CenterMask()

    def increment(self):
        return self.thisptr.Increment()
