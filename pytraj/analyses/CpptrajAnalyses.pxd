# distutils: language = c++
from pytraj.analyses.Analysis cimport _Analysis, Analysis, RetType
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core._FunctPtr cimport FunctPtr


cdef extern from "Analysis_AmdBias.h": 
    cdef cppclass _Analysis_AmdBias "Analysis_AmdBias" (_Analysis):
        _Analysis_AmdBias() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AmdBias (Analysis):
    cdef _Analysis_AmdBias* thisptr



cdef extern from "Analysis_AutoCorr.h": 
    cdef cppclass _Analysis_AutoCorr "Analysis_AutoCorr" (_Analysis):
        _Analysis_AutoCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_AutoCorr (Analysis):
    cdef _Analysis_AutoCorr* thisptr



cdef extern from "Analysis_Average.h": 
    cdef cppclass _Analysis_Average "Analysis_Average" (_Analysis):
        _Analysis_Average() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Average (Analysis):
    cdef _Analysis_Average* thisptr



cdef extern from "Analysis_Clustering.h": 
    cdef cppclass _Analysis_Clustering "Analysis_Clustering" (_Analysis):
        _Analysis_Clustering() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Clustering (Analysis):
    cdef _Analysis_Clustering* thisptr



cdef extern from "Analysis_Corr.h": 
    cdef cppclass _Analysis_Corr "Analysis_Corr" (_Analysis):
        _Analysis_Corr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Corr (Analysis):
    cdef _Analysis_Corr* thisptr



cdef extern from "Analysis_CrankShaft.h": 
    cdef cppclass _Analysis_CrankShaft "Analysis_CrankShaft" (_Analysis):
        _Analysis_CrankShaft() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrankShaft (Analysis):
    cdef _Analysis_CrankShaft* thisptr


cdef extern from "Analysis_CrdFluct.h": 
    cdef cppclass _Analysis_CrdFluct "Analysis_CrdFluct" (_Analysis):
        _Analysis_CrdFluct() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrdFluct (Analysis):
    cdef _Analysis_CrdFluct* thisptr


cdef extern from "Analysis_CrossCorr.h": 
    cdef cppclass _Analysis_CrossCorr "Analysis_CrossCorr" (_Analysis):
        _Analysis_CrossCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_CrossCorr (Analysis):
    cdef _Analysis_CrossCorr* thisptr


cdef extern from "Analysis_Divergence.h": 
    cdef cppclass _Analysis_Divergence "Analysis_Divergence" (_Analysis):
        _Analysis_Divergence() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Divergence (Analysis):
    cdef _Analysis_Divergence* thisptr



cdef extern from "Analysis_FFT.h": 
    cdef cppclass _Analysis_FFT "Analysis_FFT" (_Analysis):
        _Analysis_FFT() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_FFT (Analysis):
    cdef _Analysis_FFT* thisptr



cdef extern from "Analysis_Hist.h": 
    cdef cppclass _Analysis_Hist "Analysis_Hist" (_Analysis):
        _Analysis_Hist() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Hist (Analysis):
    cdef _Analysis_Hist* thisptr



cdef extern from "Analysis_IRED.h": 
    cdef cppclass _Analysis_IRED "Analysis_IRED" (_Analysis):
        _Analysis_IRED() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_IRED (Analysis):
    cdef _Analysis_IRED* thisptr



cdef extern from "Analysis_Integrate.h": 
    cdef cppclass _Analysis_Integrate "Analysis_Integrate" (_Analysis):
        _Analysis_Integrate() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Integrate (Analysis):
    cdef _Analysis_Integrate* thisptr



cdef extern from "Analysis_KDE.h": 
    cdef cppclass _Analysis_KDE "Analysis_KDE" (_Analysis):
        _Analysis_KDE() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_KDE (Analysis):
    cdef _Analysis_KDE* thisptr



cdef extern from "Analysis_Lifetime.h": 
    cdef cppclass _Analysis_Lifetime "Analysis_Lifetime" (_Analysis):
        _Analysis_Lifetime() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Lifetime (Analysis):
    cdef _Analysis_Lifetime* thisptr



cdef extern from "Analysis_Matrix.h": 
    cdef cppclass _Analysis_Matrix "Analysis_Matrix" (_Analysis):
        _Analysis_Matrix() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Matrix (Analysis):
    cdef _Analysis_Matrix* thisptr



cdef extern from "Analysis_MeltCurve.h": 
    cdef cppclass _Analysis_MeltCurve "Analysis_MeltCurve" (_Analysis):
        _Analysis_MeltCurve() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_MeltCurve (Analysis):
    cdef _Analysis_MeltCurve* thisptr



cdef extern from "Analysis_Modes.h": 
    cdef cppclass _Analysis_Modes "Analysis_Modes" (_Analysis):
        _Analysis_Modes() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Modes (Analysis):
    cdef _Analysis_Modes* thisptr



cdef extern from "Analysis_MultiHist.h": 
    cdef cppclass _Analysis_MultiHist "Analysis_MultiHist" (_Analysis):
        _Analysis_MultiHist() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_MultiHist (Analysis):
    cdef _Analysis_MultiHist* thisptr



cdef extern from "Analysis_Overlap.h": 
    cdef cppclass _Analysis_Overlap "Analysis_Overlap" (_Analysis):
        _Analysis_Overlap() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Overlap (Analysis):
    cdef _Analysis_Overlap* thisptr



cdef extern from "Analysis_Regression.h": 
    cdef cppclass _Analysis_Regression "Analysis_Regression" (_Analysis):
        _Analysis_Regression() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Regression (Analysis):
    cdef _Analysis_Regression* thisptr



cdef extern from "Analysis_RemLog.h": 
    cdef cppclass _Analysis_RemLog "Analysis_RemLog" (_Analysis):
        _Analysis_RemLog() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RemLog (Analysis):
    cdef _Analysis_RemLog* thisptr



cdef extern from "Analysis_Rms2d.h": 
    cdef cppclass _Analysis_Rms2d "Analysis_Rms2d" (_Analysis):
        _Analysis_Rms2d() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Rms2d (Analysis):
    cdef _Analysis_Rms2d* thisptr



cdef extern from "Analysis_RmsAvgCorr.h": 
    cdef cppclass _Analysis_RmsAvgCorr "Analysis_RmsAvgCorr" (_Analysis):
        _Analysis_RmsAvgCorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RmsAvgCorr (Analysis):
    cdef _Analysis_RmsAvgCorr* thisptr



cdef extern from "Analysis_RunningAvg.h": 
    cdef cppclass _Analysis_RunningAvg "Analysis_RunningAvg" (_Analysis):
        _Analysis_RunningAvg() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_RunningAvg (Analysis):
    cdef _Analysis_RunningAvg* thisptr


cdef extern from "Analysis_Spline.h": 
    cdef cppclass _Analysis_Spline "Analysis_Spline" (_Analysis):
        _Analysis_Spline() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Spline (Analysis):
    cdef _Analysis_Spline* thisptr



cdef extern from "Analysis_Statistics.h": 
    cdef cppclass _Analysis_Statistics "Analysis_Statistics" (_Analysis):
        _Analysis_Statistics() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Statistics (Analysis):
    cdef _Analysis_Statistics* thisptr


cdef extern from "Analysis_Timecorr.h": 
    cdef cppclass _Analysis_Timecorr "Analysis_Timecorr" (_Analysis):
        _Analysis_Timecorr() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_Timecorr (Analysis):
    cdef _Analysis_Timecorr* thisptr


cdef extern from "Analysis_VectorMath.h": 
    cdef cppclass _Analysis_VectorMath "Analysis_VectorMath" (_Analysis):
        _Analysis_VectorMath() 
        _DispatchObject * Alloc() 
        void Help()


cdef class Analysis_VectorMath (Analysis):
    cdef _Analysis_VectorMath* thisptr
