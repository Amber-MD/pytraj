# distutils: language = c++
from libcpp.string cimport string
from pytraj.trajs.Trajin cimport _Trajin, Trajin
from pytraj._FArray cimport _FArray, _FArray_iter
from pytraj.FrameArray cimport FrameArray
from pytraj.Frame cimport Frame
from pytraj.trajs.TrajectoryFile cimport _TrajectoryFile

cdef extern from "Trajin_Multi.h": 
    # Trajin_Multi.h
    ctypedef enum TargetType "Trajin_Multi::TargetType":
        pass
    cdef cppclass _Trajin_Multi "Trajin_Multi" (_Trajin):
        _Trajin_Multi() 
        void EnsembleInfo() const 
        int EnsembleSetup(_FArray&)
        int GetNextEnsemble(_FArray&)
        int EnsembleSize() const 
        # we don't need MPI here
        int EnsemblePosition(int member) const 
        bint BadEnsemble() const 
        TargetType TargetMode() const 
        string FinalCrdIndices() const 

cdef class Trajin_Multi (Trajin):
    cdef _Trajin_Multi* thisptr
    # hold iternal FrameArray
    cdef _FArray _farray
