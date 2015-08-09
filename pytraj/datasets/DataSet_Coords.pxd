# distutils: language = c++
from ..Topology cimport _Topology, Topology
from .DataSet_1D cimport _DataSet_1D, DataSet_1D
from .DataSet cimport _DataSet, DataSet, DataType
from pytraj.Frame cimport _Frame, Frame
from pytraj.AtomMask cimport _AtomMask, AtomMask


cdef extern from "DataSet_Coords.h": 
    cdef cppclass _DataSet_Coords "DataSet_Coords" (_DataSet):
        _DataSet_Coords() 
        _DataSet_Coords(DataType)
        #virtual ~_DataSet_Coords() 
        _Frame AllocateFrame() const 
        
        # virtual methods
        void AddFrame(const _Frame&) 
        void SetCRD(int, const _Frame&) 
        void GetFrame(int, _Frame&) 
        void GetFrame(int, _Frame&, const _AtomMask&) 
        # end virtual methods

        void SetTopology(const _Topology&)
        inline const _Topology& Top() const 


cdef class DataSet_Coords (DataSet):
    # DataSet has baseptr0
    cdef _DataSet_Coords* baseptr_1
    cdef Topology _top
    cdef bint py_free_mem

    # use tmpfarray object to hold Frame or Trajectory 
    # (if we want to use dset[0][0] correctly)
    cdef object tmpfarray
