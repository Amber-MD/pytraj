# distuils: language = c++

cdef class FunctPtr:
    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass

    #def get_action(self):
    #    cdef Action act = Action()
    #    act.thisptr = <_Action*> self.ptr()
    #    return act
