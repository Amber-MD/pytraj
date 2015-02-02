# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from BufferedLine cimport *


cdef extern from "Action_NMRrst.h": 
#    cdef cppclass _Action_NMRrst::Site "Action_NMRrst::Site" (_Action):
#        Site() 
#        Site(int r, const Iarray& i)
#        int ResNum() const 
#        unsigned int Nindices() const 
#        int Idx(unsigned int i) const 
#        int Count(unsigned int i) const 
#        void Increment(int c)
#        string SiteLegend(const _Topology&) const 


    cdef cppclass _Action_NMRrst "Action_NMRrst" (_Action):
        _Action_NMRrst() 
        _DispatchObject * Alloc() 
        void Help() 


#    cdef cppclass _Action_NMRrst::NOEtype "Action_NMRrst::NOEtype" (_Action):
#        NOEtype() 
#        NOEtype(const Site& s1, const Site& s2, _DataSet * d, const string& l)
#        const Site& Site1() const 
#        const Site& Site2() const 
#        double R6_Avg() const 
#        const char * legend() const 
#        _DataSet * Data() 
#        void SetR6Avg(double r6)
#        string PrintNOE() const 
#        void UpdateNOE(int i, double d2, unsigned int c1, unsigned int c2)
#        inline bint operator[( const NOEtype& rhs) const 


#cdef class Action_NMRrst::noeDataType (Action):
#    cdef _Action_NMRrst::noeDataType* thisptr
#
#cdef class Action_NMRrst::Site (Action):
#    cdef _Action_NMRrst::Site* thisptr

cdef class Action_NMRrst (Action):
    cdef _Action_NMRrst* thisptr

#cdef class Action_NMRrst::NOEtype (Action):
#    cdef _Action_NMRrst::NOEtype* thisptr

