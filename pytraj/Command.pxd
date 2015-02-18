# distutils: language = c++
from pycpptraj.CpptrajState cimport *
from pycpptraj.ArgList cimport *
from pycpptraj.DispatchObject cimport *
from pycpptraj._FunctPtr cimport FunctPtr


cdef extern from "Command.h": 
    ctypedef enum RetType "Command::RetType":
        C_OK "Command::C_OK"
        C_ERR "Command::C_ERR"
        C_QUIT "Command::C_QUIT"
    ctypedef enum CommandType "Command::CommandType":
        NONE "Command::NONE"
        PARM "Command::PARM"
        TRAJ "Command::TRAJ"
        ACTION "Command::ACTION"
        ANALYSIS "Command::ANALYSIS"
        GENERAL "Command::GENERAL"
        DEPRECATED "Command::DEPRECATED"
    ctypedef DispatchAllocatorType AllocType
    #ctypedef RetType (*CommandFxnType)(_CpptrajState&, _ArgList&, _AllocType)
    #ctypedef (*CommandHelpType)()
    #ctypedef const char* CommandKeywordType
    ctypedef struct Token:
        pass
        #CommandType Type
        #CommandKeywordType Cmd
        #AllocType Alloc
        #CommandHelpType Help
        #CommandFxnType Fxn
    ctypedef Token* TokenPtr
    cdef cppclass _Command "Command":
        void ListCommands(CommandType)
        TokenPtr SearchTokenType(_CommandType, const _ArgList& argIn)
        TokenPtr SearchToken(_ArgList&)
        RetType Dispatch(_CpptrajState&, const string&)
        RetType ProcessInput(_CpptrajState&, const string&)
        const Token& CmdToken(int idx)

cdef class Command:
    cdef _Command* thisptr

#  this "LoadCrd" method is in "Command.cpp" in cpptraj but not in the header file
# How can I include here?
#cdef extern from *:
#    #load_crd "LoadCrd"(_CpptrajState& State, _ArgList& argIn, AllocType Alloc)
#    LoadCrd(_CpptrajState& State, _ArgList& argIn, AllocType Alloc)
