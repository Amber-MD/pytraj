# distutils: language = c++
from ComplexArray cimport *


cdef extern from "PubFFT.h": 
    cdef cppclass _PubFFT "PubFFT":
        _PubFFT() 
        #~_PubFFT() 
        _PubFFT(int)
        _PubFFT(const _PubFFT&)
        #_PubFFT& operator =(const _PubFFT&)
        int size() const 
        void Forward(_ComplexArray&)
        void Back(_ComplexArray&)
        int SetupFFT_NextPowerOf2(int)
        int SetupFFTforN(int)
