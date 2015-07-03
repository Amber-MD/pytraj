# distutils: language = c++
from __future__ import absolute_import
from .CpptrajState cimport *
from ..ArgList cimport *
from .DispatchObject cimport *
from ._FunctPtr cimport FunctPtr


cdef extern from "Command.h": 
    ctypedef enum RetType "Command::RetType":
        pass
    cdef cppclass _Command "Command":
        @staticmethod
        RetType ProcessInput(_CpptrajState&, const string&)

cdef class Command:
    cdef _Command* thisptr

    def __cinit__(self):
        self.thisptr = new _Command()

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def get_state(cls, trajin_text):
        cdef CpptrajState cppstate = CpptrajState()
        trajin_text = trajin_text.encode()
        _Command.ProcessInput(cppstate.thisptr[0], trajin_text)
        return cppstate
