# distutils: language = c++


cdef class Command:
    def __cinit__(self):
        self.thisptr = new _Command()

    def __dealloc__(self):
        del self.thisptr

    #def ListCommands(self, cmtype=None):
    #    self.thisptr.ListCommands(self.mode('CommandType')[cmtype])

    #def TokenPtr SearchTokenType(self,CommandType cmmode, ArgList argIn):
    #def TokenPtr SearchToken(self,ArgList):

    #def Dispatch(self, CpptrajState state, string cmdin):
    #    return self.thisptr.Dispatch(state.thisptr[0], cmdin)

    def process_input(self, CpptrajState cppstate, string fnameIn):
        return self.thisptr.ProcessInput(cppstate.thisptr[0], fnameIn)
            
    #cdef Token&& CmdToken(self, int idx):
    #    return self.thisptr.CmdToken(idx)

    #@classmethod
    #def load_crd(cls, CpptrajState state, ArgList arglist, obj):
    #    """temp doc: load_crd(self, CpptrajState state, ArgList arglist, obj)
    #    obj :: action or analysis object
    #    """
    #    cdef FunctPtr func = <FunctPtr> obj.alloc()
    #    _Command.LoadCrd(state.thisptr[0], arglist.thisptr[0], func.ptr)
