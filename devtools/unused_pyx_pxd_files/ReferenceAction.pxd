# distutils: language = c++
from __future__ import absolute_import
from libcpp.string cimport string
from pytraj.trajs.Trajin_Single cimport *
from pytraj.Frame cimport _Frame, Frame
from pytraj.Vec3 cimport _Vec3, Vec3
from pytraj.Topology cimport _Topology, Topology
from ..ReferenceFrame cimport *


cdef extern from "ReferenceAction.h": 
    cdef cppclass _ReferenceAction "ReferenceAction":
        _Reference_Action() 
        void SetRefStructure(const _Frame&, bint, bint)
        int InitRef(bint, bint, bint, bint, const string&, const _Reference_Frame&,
                    _Topology *, const string&, _ArgList&, const char *)
        int SetupRef(const _Topology&, int, const char *)
        void _ActionRef(const _Frame&, bint, bint)
        bint Previous()
        const char * RefModeString()
        const _Frame& Ref_Frame()
        const _Frame& SelectedRef()
        const _Vec3& RefTrans()


cdef class ReferenceAction:
    cdef _ReferenceAction* thisptr
