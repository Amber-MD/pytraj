from cython cimport view

def get_positive_idx(idx, size):
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
