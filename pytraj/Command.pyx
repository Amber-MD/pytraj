# distutils: language = c++


cdef class Command:
    def __cinit__(self):
        self.thisptr = new _Command()

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def process_input(cls, CpptrajState cppstate, fnameIn):
        fnameIn = fnameIn.encode()
        _Command.ProcessInput(cppstate.thisptr[0], fnameIn)
