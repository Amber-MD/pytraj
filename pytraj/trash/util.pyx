from libcpp.vector cimport vector
from ParameterTypes cimport ptype

cdef vector[ptype*] convert_list_to_vector(list listin):
    cdef vector[ptype*] v
    cdef ptype* ptr

    for ins in listin:
        ptr = <ptype*>ins.thisptr
        v.push_back(ptr)
    return v

