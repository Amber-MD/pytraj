# distutils: language = c++
from libcpp.string cimport string
from ..core.cpptraj_core cimport (_DispatchObject, DispatchObject,  FunctPtr)
from ..datafiles.datafiles cimport  _DataFileList, DataFileList
from ..core.TopologyList cimport _TopologyList, TopologyList
from ..core.cpptraj_core cimport _ArgList, ArgList
from ..datasets.DataSetList cimport _DataSetList, DataSetList
from ..Topology cimport _Topology, Topology
from ..Frame cimport _Frame, Frame


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
