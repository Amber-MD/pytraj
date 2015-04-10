# distutils: language = c++


cdef class Command:
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
