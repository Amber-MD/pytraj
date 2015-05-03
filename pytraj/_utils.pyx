def _list2d_to_memview(mylist):
    """convert 2D list of integers to cython.view.array
    """
    cdef int n_counts = len(mylist)
    cdef int n_items = len(mylist[0])
    cdef view.array arr0 = view.array(shape=(n_counts, n_items), 
                           itemsize=sizeof(int), format="i")
    cdef int i, j
    for i in range(arr0.shape[0]):
        for j in range(arr0.shape[1]):
            arr0[i, j] = mylist[i][j]
    return arr0
