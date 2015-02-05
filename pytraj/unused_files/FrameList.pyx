# distutils: language = c++


cdef class FrameList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _FrameList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    # TODO : do we really need this?

    #def clear(self):
    #    self.thisptr.Clear()

    #def set_debug(self, int idx):
    #    """
    #    Input:
    #    ====
    #    idx : int
    #    """
    #    self.thisptr.SetDebug(idx)

    #def get_active_reference(self):
    #    """
    #    return Frame instance
    #    """
    #    cdef Frame frame = Frame()
    #    frame.thisptr[0] = self.thisptr.ActiveReference()
    #    return frame

    #def set_active_ref(self, int numIn):
    #    """
    #    set_active_ref for frame with idx 

    #    Input:
    #    =====
    #    idx :: int
    #    """
    #    return self.thisptr.SetActiveRef(numIn)

    #def add_reference(self, ArgList arglist, TopologyList toplist):
    #    """
    #    Input:
    #    =====
    #    arglist :: ArgList instance
    #    toplist :: TopologyList instance
    #    """
    #    return self.thisptr.AddRefFrame(arglist.thisptr[0], toplist.thisptr[0])

    #def get_frame_from_args(self,ArgList arglist):
    #    """
    #    Input:
    #    =====
    #    arglist :: ArgList instance:w
    #    """
    #    cdef ReferenceFrame ref = ReferenceFrame()
    #    ref.thisptr[0] = self.thisptr.GetFrameFromArgs(arglist.thisptr[0])
    #    return ref

    #def get_frame_by_name(self, string name):
    #    """
    #    Input:
    #    =====
    #    name :: str
    #    """
    #    cdef ReferenceFrame ref = ReferenceFrame()
    #    ref.thisptr[0] = self.thisptr.GetFrameByName(name)
    #    return ref

    #def list(self):
    #    self.thisptr.List()

    #@property
    #def n_frames(self):
    #    return self.thisptr.NumFrames()
