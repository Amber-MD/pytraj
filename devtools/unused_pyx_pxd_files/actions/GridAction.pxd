# distutils: language = c++
from pytraj.datasets.DataSetList cimport _DataSetList, DataSetList
from pytraj.datasets.DatasetGridFloat cimport _DatasetGridFloat, DatasetGridFloat
from pytraj.Frame cimport _Frame, Frame
from pytraj.ArgList cimport _ArgList, ArgList
from pytraj.AtomMask cimport _AtomMask, AtomMask
from pytraj.Topology cimport _Topology, Topology


cdef extern from "GridAction.h": 
    # GridAction.h
    ctypedef enum GridModeType "GridAction::GridModeType":
        ORIGIN "GridAction::ORIGIN"
        BOX "GridAction::BOX"
        MASKCENTER "GridAction::MASKCENTER"
        SPECIFIEDCENTER "GridAction::SPECIFIEDCENTER"
    cdef cppclass _GridAction "GridAction":
        Grid_Action() 
        _DatasetGridFloat * GridInit(const char *, _ArgList&, _DataSetList&)
        void GridInfo(const _DatasetGridFloat&)
        int GridSetup(const _Topology&)
        inline void GridFrame(const _Frame&, const _AtomMask&, _DatasetGridFloat&)
        GridModeType GridMode() const 
        const _AtomMask& CenterMask() const 
        float Increment() const 


cdef class GridAction:
    cdef _GridAction* thisptr
