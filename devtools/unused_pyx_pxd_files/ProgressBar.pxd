# distutils: language = c++


cdef extern from "ProgressBar.h": 
    cdef cppclass _ParallelProgress "ParallelProgress":
        _ParallelProgress() 
        _ParallelProgress(int m)
        _ParallelProgress(const _ParallelProgress &)
        void SetThread(int t)
        void Update(int it)
        void Finish() 
