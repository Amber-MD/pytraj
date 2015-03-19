from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as incr
from pytraj.Frame cimport Frame
from pytraj.CpptrajStdio cimport SetWorldSilent as cpptraj_SetWorldSilent

def set_world_silent(turnoff=True):
    cpptraj_SetWorldSilent(turnoff)

def _tease_FrameArray(farray, traj, int nloop):
    cdef Frame frame
    cdef int i

    for i in range(nloop):
        for frame in traj:
            farray.append(frame, copy=False)

# STATUS: not worked yet
# wrong memory (got all 0.0 )
cdef double[:] vec_to_memview_double(vector[double] varray):
    cdef int idx = varray.size()
    cdef double* double_ptr
    cdef cyarray arr0

    double_ptr = &varray[0]
    arr0 = <double[:idx]> double_ptr
    return arr0

# STATUS: not worked yet
# wrong memory (got all 0.0 )
cdef int[:] vec_to_memview_int(vector[int] varray):
    cdef int idx = varray.size()
    cdef int* int_ptr
    cdef cyarray arr0

    int_ptr = &varray[0]
    arr0 = <int[:idx]> int_ptr
    return arr0

def test():
    cdef vector[double] v = range(100)
    print ("vector v = ", v)
    return vec_to_memview_double(v)
