# distutils: language = c++
from pytraj.CpptrajState cimport *
from pytraj.ArgList cimport *
from pytraj.DispatchObject cimport *
from pytraj._FunctPtr cimport FunctPtr


cdef extern from "src/Command.h": 
    ctypedef enum RetType "Command::RetType":
        pass
    cdef cppclass _Command "Command":
        @staticmethod
        RetType ProcessInput(_CpptrajState&, const string&)

cdef class Command:
    cdef _Command* thisptr
