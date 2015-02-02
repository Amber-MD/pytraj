# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject
#from ImagedAction cimport *
#from DataSet_MatrixDbl cimport *
#from DataSet_integer cimport *


cdef extern from "Action_NativeContacts.h": 
    cdef cppclass _Action_NativeContacts "Action_NativeContacts" (_Action):
        _Action_NativeContacts() 
        _DispatchObject * Alloc() 
        void Help() 


#    cdef cppclass _Action_NativeContacts::contactType "Action_NativeContacts::contactType" (_Action):
#        contactType() 
#        contactType(const string& id)
#        const char * id() const 
#        int Nframes() const 
#        double Avg() const 
#        double Stdev() const 
#        _DataSet_integer& Data() 
#        void Increment(int fnum, double d, double d2)
#        void Finalize() 
#        bint operator[( const contactType& rhs) const 
#        void SetData(_DataSet * ds)


cdef class Action_NativeContacts (Action):
    cdef _Action_NativeContacts* thisptr
