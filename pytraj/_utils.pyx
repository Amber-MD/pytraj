from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as incr
from pytraj.Frame cimport Frame

def _tease_FrameArray(farray, traj, int nloop):
    cdef Frame frame
    cdef int i

    for i in range(nloop):
        for frame in traj:
            farray.append(frame, copy=False)
