# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from FusedType cimport *

cdef generator_template(self, CommonClassCpptraj newins, IteratorCpptraj it, begin_it, end_it):

    it = begin_it()
    while it != end_it():
        newins  = class_name()
        newins.thisptr[0] = deref(it)
        yield newins
        incr(it)
        

