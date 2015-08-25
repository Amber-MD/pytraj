
from .Frame cimport Frame, _Frame
from .core.Box cimport _Box, Box

def _fast_iterptr(double[:, :, :] xyz, int n_atoms):
    '''noxbox
    '''
    cdef int i
    cdef int n_frames = xyz.shape[0]

    # just create a pointer
    cdef Frame frame = Frame(n_atoms, xyz[0], _as_ptr=True)
    cdef _Frame _frame

    for i in range(n_frames):
        frame.thisptr.SetXptr(n_atoms, &xyz[i, 0, 0])
        yield frame

def _fast_iterptr_withbox(double[:, :, :] xyz, double[:, :] boxes, int n_atoms):
    # withbox
    cdef int i
    cdef int n_frames = xyz.shape[0]

    # just create a pointer
    cdef Frame frame = Frame(n_atoms, xyz[0], _as_ptr=True)
    cdef _Frame _frame

    for i in range(n_frames):
        frame.thisptr.SetXptr(n_atoms, &xyz[i, 0, 0])
        frame.thisptr.SetBox(_Box(&boxes[i, 0]))
        yield frame
