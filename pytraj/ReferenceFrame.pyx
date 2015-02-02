# distutils: language = c++


cdef class ReferenceFrame:
    """
    Original cpptraj doc:
    ====================
        Hold single frame coordinates to be used as a reference.
        The frame and topology are stored as pointers instead of classes to
        save memory; this way multiple actions can use the same reference 
        structure without each having to have a different copy. Because
        of this, memory is not freed in ReferenceFrame destructor to avoid
        potential double-frees when ReferenceFrame is used in e.g. vectors.
        Freeing must be accomplished with the ClearRef function.
    """
    def __cinit__(self):
        self.thisptr = new _ReferenceFrame()

    def __dealloc__(self):
        del self.thisptr

    # TODO : do we really need this module?

    #@property
    #def frame(self):
    #    """return Frame instance"""
    #    cdef Frame frame = Frame()
    #    frame.thisptr[0] = self.thisptr.Coord()
    #    return frame

    #@property
    #def top(self):
    #    cdef Topology top = Topology()
    #    top.thisptr[0] = self.thisptr.Parm()
    #    return top

    #def is_empty(self):
    #    return self.thisptr.empty()

    #def frame_name(self):
    #    """Return FileName instance
    #    TODO: Should we really need FileName class?
    #    """
    #    cdef FileName fn = FileName()
    #    fn.thisptr[0] = self.thisptr.FrameName()
    #    return fn

    #def load_ref(self, string filename="", Topology top=Topology(), int debug=0, ArgList arglist=ArgList(), string mask=""):
    #    """Temp doc: load_ref(self, string filename, Topology parmIn, debug=0, *args)"""
    #    return self.thisptr.LoadRef(filename, arglist.thisptr[0], top.thisptr, mask, debug)

    #def strip_ref(self, AtomMask atm):
    #    """temp doc: strip_ref(self, AtomMask atm)"""
    #    return self.thisptr.StripRef(atm.thisptr[0])

    #def info(self):
    #    self.thisptr.RefInfo()

    #def clear_ref(self):
    #    "Free memory"
    #    self.thisptr.ClearRef()
