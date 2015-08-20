
from .Frame cimport Frame, _Frame

def _fastiter(self, int n_atoms):
    cdef double[:, :, ::1] view = self.xyz
    cdef int i
    cdef int n_frames = view.shape[0]

    # just create a pointer
    cdef Frame frame = Frame()
    cdef _Frame _frame

    #frame._as_view = True

    for i in range(n_frames):
        #_frame = _Frame(n_atoms, &view[i, 0, 0])
        frame = Frame(n_atoms, view[i], _as_ptr=True)
        yield frame
        #del _frame # don't deallocate, will get double-free

def _fastiter_ptr(self, int n_atoms):
    cdef double[:, :, ::1] view = self.xyz
    cdef int i
    cdef int n_frames = view.shape[0]

    # just create a pointer
    cdef Frame frame = Frame(n_atoms, view[0], _as_ptr=True)
    cdef _Frame _frame

    #frame._as_view = True

    for i in range(n_frames):
        #_frame = _Frame(n_atoms, &view[i, 0, 0])
        frame.thisptr.SetXptr(&view[i, 0, 0])
        yield frame
        #del _frame # don't deallocate, will get double-free
