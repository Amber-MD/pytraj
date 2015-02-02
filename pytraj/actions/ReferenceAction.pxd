# distutils: language = c++
from libcpp.string cimport string
from Trajin_Single cimport *
from ReferenceFrame cimport *


cdef extern from "ReferenceAction.h": 
    cdef cppclass _ReferenceAction "ReferenceAction":
        _Reference_Action() 
        inline void SetRefStructure(const _Frame&, bint, bint)
        int InitRef(bint, bint, bint, bint, const string&, _ReferenceFrame&, _Topology *, const string&, _ArgList&, const char *)
        int SetupRef(const _Topology&, int, const char *)
        inline void _ActionRef(const _Frame&, bint, bint)
        bint Previous() const 
        const char * RefModeString() const 
        const _Frame& RefFrame() const 
        const _Frame& SelectedRef() const 
        const _Vec3& RefTrans() const 
