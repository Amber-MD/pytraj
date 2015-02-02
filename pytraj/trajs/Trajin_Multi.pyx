# distutils: language = c++
from pytraj.cpptraj_dict import get_key, TargetDict

cdef class Trajin_Multi (Trajin):
    def __cinit__(self):
        self.thisptr = new _Trajin_Multi()
        self.baseptr1 = <_Trajin*> self.thisptr
        self.baseptr_0 = <_TrajectoryFile*> self.thisptr

    def __dealloc__(self):
        del self.thisptr

    def ensemble_info(self):
        self.thisptr.EnsembleInfo()

    def ensemble_setup(self, FrameArray2 farray2):
        return self.thisptr.EnsembleSetup(farray2.thisptr[0])

    def get_next_ensemble(self, FrameArray2 farray2):
        return self.thisptr.GetNextEnsemble(farray2.thisptr[0])

    @property
    def size(self):
        return self.thisptr.EnsembleSize()

    def ensemble_frame_num(self):
        return self.thisptr.EnsembleFrameNum()

    def ensemble_position(self, int memberidx):
        return self.thisptr.EnsemblePosition(memberidx)

    def bad_ensemble(self):
        return self.thisptr.BadEnsemble()

    #def TargetType TargetMode(self):
    def target_mode(self):
        return get_key(self.thisptr.TargetMode(), TargetDict)

    def final_crd_indices(self):
        return self.thisptr.FinalCrdIndices()

