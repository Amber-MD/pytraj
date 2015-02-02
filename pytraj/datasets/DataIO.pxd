# distutils: language = c++
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.datasets.DataSet cimport _DataSet, DataSet
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.CpptrajFile cimport _CpptrajFile, CpptrajFile
from pytraj.BaseIOtype cimport _BaseIOtype, BaseIOtype


cdef extern from "DataIO.h": 
    #abstract class
    cdef cppclass _DataIO "DataIO":
        _DataIO()
        _DataIO(bint v1, bint v2, bint v3)
        bint CheckValidFor(const _DataSet&)const 
