cimport cython
from .Frame cimport Frame

@cython.boundscheck(False)
@cython.wraparound(False)
def _faster_iter(double[:, :, :] xyz, int n_atoms, int n_frames):
    cdef Frame frame
    cdef int i
    cdef double* ptr
    cdef double[:, :] aview

    for i in range(n_frames):
        frame = Frame(n_atoms)
        ptr = frame.thisptr.xAddress()
        aview = <double[:n_atoms, :3]> ptr
        aview[:] = xyz[i]
        yield frame
