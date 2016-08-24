cimport cython
from cython cimport view
from .frame cimport Frame, _Frame
from .core.box cimport _Box, Box
from .topology cimport Topology


__all__ = ['_fast_count', 'get_positive_idx', '_int_array1d_like_to_memview',
           '_int_array2d_like_to_memview',
           '_fast_iterptr', '_fast_iterptr_withbox', ]


def _fast_count(cython.integral[:] values, int target):
    cdef int i
    cdef int count = 0

    for i in range(values.shape[0]):
        if values[i] == target:
            count += 1
    return count


def get_positive_idx(int idx, int size):
    """Used for negative indexing"""
    if idx < 0:
        idx = size + idx
        if idx < 0:
            raise ValueError("index is out of range")
    if idx >= size:
        raise ValueError("index is out of range")
    return idx


def _int_array2d_like_to_memview(mylist):
    """convert 2D array-like of integers to cython.view.array
    """
    cdef int n_counts = len(mylist)
    cdef int n_items = len(mylist[0])
    cdef int[:, :] int_view
    cdef view.array arr0 = view.array(shape=(n_counts, n_items),
                                      itemsize=sizeof(int), format="i")
    cdef int i, j
    try:
        int_view = mylist
    except:
        for i in range(arr0.shape[0]):
            for j in range(arr0.shape[1]):
                arr0[i, j] = mylist[i][j]
        int_view = arr0
    return int_view


def _int_array1d_like_to_memview(mylist):
    """convert 1D array-like of integers to cython.view.array
    """
    cdef int n_counts = len(mylist)
    cdef int[:] int_view
    cdef view.array arr0 = view.array(shape=(n_counts,),
                                      itemsize=sizeof(int), format="i")
    cdef int i

    try:
        # if `mylist` has buffer
        int_view = mylist
    except:
        # if not, try to make a copy
        for i in range(arr0.shape[0]):
            arr0[i] = mylist[i]
        int_view = arr0
    return int_view


def _fast_iterptr(double[:, :, :] xyz, int n_atoms, indices, Topology topology=None,
                  double[:] mass=None):
    '''noxbox
    '''
    cdef int i
    cdef int n_frames = xyz.shape[0]

    # just create a pointer
    cdef Frame frame = Frame(n_atoms, xyz[0], _as_ptr=True)
    if topology is not None:
        frame.set_mass(topology)
    elif mass is not None:
        frame._set_mass_from_array(mass)

    for i in indices:
        frame.thisptr.SetXptr(n_atoms, & xyz[i, 0, 0])
        yield frame


def _fast_iterptr_withbox(double[:, :, :] xyz, double[:, :] boxes, int n_atoms, indices,
                          Topology topology=None, double[:] mass=None):
    # withbox
    cdef int i
    cdef int n_frames = xyz.shape[0]

    # just create a pointer
    cdef Frame frame = Frame(n_atoms, xyz[0], _as_ptr=True)
    if topology is not None:
        frame.set_mass(topology)
    elif mass is not None:
        frame._set_mass_from_array(mass)

    for i in indices:
        frame.thisptr.SetXptr(n_atoms, & xyz[i, 0, 0])
        frame.thisptr.SetBox(_Box( & boxes[i, 0]))
        yield frame
