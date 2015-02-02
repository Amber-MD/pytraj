# distutils: language = c++


cdef extern from "Random.h": 
    cdef cppclass _Random_Number "Random_Number":
        _Random_Number() 
        void rn_set(int)
        void rn_set() 
        double rn_gen() 
        double rn_gauss(double, double)
        bint IsSet() const 
