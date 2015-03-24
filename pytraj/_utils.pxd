# distutils: language = c++
from libc.stdlib cimport malloc 
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython cimport view
cimport cython
from cython.view cimport array as cyarray
from cpython.array import array as pyarray

num_fused_type = cython.fused_type(int, float, double)

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

# should use fused_type: TODO: how?
# STATUS: not worked yet. 
# wrong memory (got all 0.0 )
#cdef inline num_fused_type[:] vec_to_memview(vec_to_memview[num_fused_type] varray, string dtype=?):
#    cdef int idx = varray.size()
#    cdef int* int_ptr
#    cdef double* d_ptr
#
#    if dtype == 'int':
#        int_ptr = &varray[0]
#        return <int[:varray.size()]> int_ptr
#    elif dtype == 'double':
#        d_ptr = &varray[0]
#        return <double[:varray.size()]> d_ptr
#    else:
#        raise ValueError("only support casting for int or double")

# wrong memory (got all 0.0 )
cdef inline double[:] vec_to_memview_double(vector[double] varray)
    #cdef int idx = varray.size()
    #cdef double* double_ptr
    #cdef cyarray arr0

    #double_ptr = &varray[0]
    #arr0 = <double[:idx]> double_ptr
    #return arr0

# wrong memory (got all 0.0 )
cdef inline int[:] vec_to_memview_int(vector[int] varray)
    #cdef int idx = varray.size()
    #cdef int* int_ptr
    #cdef cyarray arr0

    #int_ptr = &varray[0]
    #arr0 = <int[:idx]> int_ptr
    #return arr0
