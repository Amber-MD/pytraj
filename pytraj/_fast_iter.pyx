
from .Frame cimport Frame, _Frame

def _fastiter(self, int n_atoms):
    cdef double[:, :, ::1] view = self.xyz
    cdef int i
    cdef int n_frames = view.shape[0]

    # just create a pointer
    cdef Frame frame = Frame()
    cdef _Frame* _frame

    for i in range(n_frames):
        _frame = new _Frame(n_atoms, &view[i, 0, 0])
        frame.thisptr = _frame
        yield frame
        #del _frame # don't deallocate, will get double-free
