# distutils: language = c++
from libcpp.string cimport string
from pytraj.core.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.core.DataFileList cimport _DataFileList, DataFileList
from pytraj.core.TopologyList cimport _TopologyList, TopologyList
from pytraj.core._FunctPtr cimport FunctPtr
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.datasets.DataSetList cimport _DataSetList, DataSetList
from pytraj.Topology cimport _Topology, Topology
from pytraj.Frame cimport _Frame, Frame


cdef extern from "Analysis.h": 
    # Analysis.h
    ctypedef enum RetType "Analysis::RetType":
        pass
    cdef cppclass _Analysis "Analysis":
        #virtual ~_Analysis() 
        RetType Setup(_ArgList&, _DataSetList *, _TopologyList *, _DataFileList *, int)
        RetType Analyze()


cdef class Analysis:
    cdef _Analysis* baseptr
