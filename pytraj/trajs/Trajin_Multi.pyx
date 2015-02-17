# distutils: language = c++

from cython.operator cimport dereference as deref, preincrement as incr
from pytraj.cpptraj_dict import get_key, TargetDict

cdef class Trajin_Multi (Trajin):
    def __cinit__(self):
        self.thisptr = new _Trajin_Multi()
        self.baseptr_1 = <_Trajin*> self.thisptr
        self.baseptr0 = <_TrajectoryFile*> self.thisptr

    def __dealloc__(self):
        del self.thisptr

    def ensemble_info(self):
        self.thisptr.EnsembleInfo()

    def ensemble_setup(self):
        self.thisptr.EnsembleSetup(self._farray)

    def get_next_ensemble(self, FrameArray farray):
        cdef _FArray_iter it 
        cdef Frame frame = Frame()

        self.thisptr.GetNextEnsemble(self._farray)
        it = self._farray.begin()
        while it != self._farray.end():
            frame.thisptr = &(deref(it))
            farray.append(frame, copy=False)
            incr(it)

    @property
    def n_ensembles(self):
        return self.thisptr.EnsembleSize()

    def ensemble_position(self, int memberidx):
        return self.thisptr.EnsemblePosition(memberidx)

    def bad_ensemble(self):
        return self.thisptr.BadEnsemble()

    @property
    def target_mode(self):
        return get_key(self.thisptr.TargetMode(), TargetDict)

    def final_crd_indices(self):
        return self.thisptr.FinalCrdIndices().decode()
