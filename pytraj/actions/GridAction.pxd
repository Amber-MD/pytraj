# distutils: language = c++
from pytraj.DataSetList cimport _DataSetList, DataSetList
from pytraj.FrameList cimport _FrameList, FrameList
from pytraj.datasets.DataSet_GridFlt cimport _DataSet_GridFlt, DataSet_GridFlt
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
        _DataSet_GridFlt * GridInit(const char *, _ArgList&, _DataSetList&, const _FrameList&)
        void GridInfo(const _DataSet_GridFlt&)
        int GridSetup(const _Topology&)
        inline void GridFrame(const _Frame&, const _AtomMask&, _DataSet_GridFlt&)
        GridModeType GridMode() const 
        const _AtomMask& CenterMask() const 
        float Increment() const 


cdef class GridAction:
    cdef _GridAction* thisptr

