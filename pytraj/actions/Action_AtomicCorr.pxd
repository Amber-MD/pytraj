# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_AtomicCorr.h": 
    cdef cppclass _Action_AtomicCorr "Action_AtomicCorr" (_Action):
        _Action_AtomicCorr() 
        _DispatchObject * Alloc() 
        void Help() 


#    cdef cppclass _Action_AtomicCorr::AtomVector "Action_AtomicCorr::AtomVector" (_Action):
#        _AtomVector() 
#        _AtomVector(const string& sIn, int idxIn)
#        int operator -(const _AtomVector& rhs)
#        bint empty() 
#        size_t size() 
#        void push_back(float fval)
#        const string& Label() 
#        _Vec3 VXYZ(int idx)


cdef class Action_AtomicCorr (Action):
    cdef _Action_AtomicCorr* thisptr

#cdef class Action_AtomicCorr::AtomVector (Action):
#    cdef _Action_AtomicCorr::AtomVector* thisptr

