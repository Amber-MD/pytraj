# distutils: language = c++
from pytraj.ReplicaDimArray cimport _ReplicaDimArray
from pytraj.Box cimport _Box, Box


cdef extern from "CoordinateInfo.h": 
    cdef cppclass _CoordinateInfo "CoordinateInfo":
        _CoordinateInfo() 
        _CoordinateInfo(const _Box& b, bint v, bint t, bint m)
        _CoordinateInfo(const _ReplicaDimArray& r, const _Box& b, bint v, bint t, bint m, bint f)
        bint HasBox() const 
        const _Box& TrajBox() const 
        bint HasVel() const 
        bint HasTemp() const 
        bint HasTime() const 
        bint HasForce() const 
        bint HasReplicaDims() const 
        const _ReplicaDimArray& ReplicaDimensions() const 
        void SetTime(bint m)
        void SetTemperature(bint t)
        void SetVelocity(bint v)
        void SetBox(const _Box& b)

cdef class CoordinateInfo:
    cdef _CoordinateInfo* thisptr

