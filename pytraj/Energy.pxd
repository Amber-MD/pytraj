# distutils: language = c++

from libcpp.vector cimport vector
from pytraj.Topology cimport *
from pytraj.Timer cimport *
from pytraj.Frame cimport Frame, _Frame


cdef extern from "Energy.h": 
    cdef cppclass _Energy_Amber "Energy_Amber":
        _Energy_Amber() 
        double E_bond(const _Frame&, const _Topology&, const _AtomMask&)
        double E_angle(const _Frame&, const _Topology&, const _AtomMask&)
        double E_torsion(const _Frame&, const _Topology&, const _AtomMask&)
        double E_14_Nonbond(const _Frame&, const _Topology&, const _AtomMask&, double&)
        double E_Nonbond(const _Frame&, const _Topology&, const _AtomMask&, double&)
        void SetDebug(int d)
        void PrintTiming() const 


cdef class Energy_Amber:
    cdef _Energy_Amber* thisptr

