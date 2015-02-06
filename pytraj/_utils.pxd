# distutils: language = c++
from libc.stdlib cimport malloc 
from libcpp.vector cimport vector
from cython cimport view

cdef inline int get_positive_idx(idx, size):
    # TODO : do we need this method?
    # we can we memoryview to get slicing too
    """Used for negative indexing"""
    if idx < 0:
        idx = size + idx
        if idx < 0:
            raise ValueError("index is out of range")
    if idx >= size:
        raise ValueError("index is out of range")
    return idx

cdef inline char** list_to_char_pp(args):
     cdef char** c_argv
     args = [str(x) for x in args]
     c_argv = <char**>malloc(sizeof(char*) * len(args)) 
     for idx, s in enumerate(args):
         c_argv[idx] = s
     return c_argv

cdef inline _get_buffer1D(int N, double* ptr):
    cdef view.array my_arr
    my_arr = <double[:N]> ptr
    return my_arr

# from Cython website
# http://docs.cython.org/src/tutorial/strings.html
from cpython.version cimport PY_MAJOR_VERSION
cdef inline unicode _ustring(s):
    if type(s) is unicode:
        # fast path for most common case(s)
        return <unicode>s
    elif PY_MAJOR_VERSION < 3 and isinstance(s, bytes):
        # only accept byte strings in Python 2.x, not in Py3
        return (<bytes>s).decode('ascii')
    elif isinstance(s, unicode):
        # an evil cast to <unicode> might work here in some(!) cases,
        # depending on what the further processing does.  to be safe,
        # we can always create a copy instead
        return unicode(s)
    else:
        raise TypeError("")
