# distutils: language = c++
from ..ArgList cimport _ArgList, ArgList
from .DataSet cimport _DataSet, DataSet
from .DataSetList cimport _DataSetList, DataSetList
from ..core.BaseIOtype cimport _BaseIOtype, BaseIOtype


cdef extern from "DataIO.h": 
    #abstract class
    cdef cppclass _DataIO "DataIO":
        _DataIO()
        _DataIO(bint v1, bint v2, bint v3)
        bint CheckValidFor(const _DataSet&)const 
