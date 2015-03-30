# distutils: language = c++
from libcpp.string cimport string
from pytraj.analyses.Analysis cimport *
from pytraj.datasets.DataSet_2D cimport *
from pytraj.datasets.DataSet_Modes cimport *
from pytraj.DispatchObject cimport _DispatchObject, DispatchObject
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.DataFileList cimport _DataFileList, DataFileList
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.TopologyList cimport _TopologyList, TopologyList
from pytraj._FunctPtr cimport FunctPtr
from pytraj.Topology cimport _Topology, Topology
from pytraj.Frame cimport _Frame, Frame
from pytraj.FrameArray cimport FrameArray


cdef extern from "Analysis_Statistics.h": 
    cdef cppclass _Analysis_Statistics "Analysis_Statistics" (_Analysis):
        _Analysis_Statistics() 
        _DispatchObject * Alloc() 
        void Help() 
        RetType Setup(_ArgList&, _DataSetList *, _TopologyList *, _DataFileList *, int)
        RetType Analyze() 


cdef class Analysis_Statistics (Analysis):
    cdef _Analysis_Statistics* thisptr

