# distutils: language = c++

from libc.stdlib cimport malloc 
from libcpp.vector cimport vector

cdef inline char** list_to_char_pp(args):
     cdef char** c_argv
     args = [str(x) for x in args]
     c_argv = <char**>malloc(sizeof(char*) * len(args)) 
     for idx, s in enumerate(args):
         c_argv[idx] = s
     return c_argv
