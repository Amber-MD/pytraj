from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as incr
<<<<<<< HEAD
from pytraj.Frame cimport Frame
from pytraj.CpptrajStdio cimport SetWorldSilent as cpptraj_SetWorldSilent

def set_world_silent(turnoff=True):
    cpptraj_SetWorldSilent(turnoff)
=======
from pytraj.CpptrajStdio cimport SetWorldSilent as cpptraj_SetWorldSilent
>>>>>>> 018cb5b7e1cdde3dc545528846f0593eaf55bc58

def set_world_silent(turnoff=True):
    cpptraj_SetWorldSilent(turnoff)

<<<<<<< HEAD
    for i in range(nloop):
        for frame in traj:
            farray.append(frame, copy=False)

=======
>>>>>>> 018cb5b7e1cdde3dc545528846f0593eaf55bc58
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
