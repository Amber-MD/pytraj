# distutils: language = c++
from pytraj.PubFFT cimport *
from pytraj.ComplexArray cimport *


cdef extern from "Corr.h": 
    cdef cppclass _CorrF_Direct "CorrF_Direct":
        _CorrF_Direct() 
        _CorrF_Direct(int stepsIn)
        void Allocate(int)
        void AutoCorr(_ComplexArray&)
        void CrossCorr(_ComplexArray&, const _ComplexArray&)
    cdef cppclass _CorrF_FFT "CorrF_FFT":
        _CorrF_FFT() 
        _CorrF_FFT(int stepsIn)
        void Allocate(int)
        void AutoCorr(_ComplexArray&)
        void CrossCorr(_ComplexArray&, _ComplexArray&)
        _ComplexArray Array() 
