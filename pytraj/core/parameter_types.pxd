# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef vector[_BondParmType] BondParmArray
ctypedef vector[_BondType] BondArray
ctypedef vector[_AngleParmType] AngleParmArray
ctypedef vector[_AngleType] AngleArray
ctypedef vector[_DihedralParmType] DihedralParmArray
ctypedef vector[_DihedralType] DihedralArray
ctypedef vector[_NonbondType] NonbondArray


cdef extern from "ParameterTypes.h":
    # Basic type classes used in topology.pyx
    cdef cppclass _AngleType "AngleType":
        _AngleType()
        inline int A1() const
        inline int A2() const
        inline int A3() const
        inline int Idx() const

    cdef cppclass _BondType "BondType":
        _BondType()
        inline int A1() const
        inline int A2() const
        inline int Idx() const

    cdef cppclass _DihedralType "DihedralType":
        _DihedralType()
        int A1()
        int A2()
        int A3()
        int A4()
        int Idx()
        bint Skip14()
        bint IsImproper()

    # Parameter types
    cdef cppclass _BondParmType "BondParmType":
        _BondParmType()

    cdef cppclass _AngleParmType "AngleParmType":
        _AngleParmType()

    cdef cppclass _DihedralParmType "DihedralParmType":
        _DihedralParmType()

    # Other parameter types for completeness
    cdef cppclass _NonbondType "NonbondType":
        _NonbondType()

    cdef cppclass _NonbondParmType "NonbondParmType":
        _NonbondParmType()

    cdef cppclass _CapParmType "CapParmType":
        _CapParmType()

    cdef cppclass _LES_ParmType "LES_ParmType":
        _LES_ParmType()

    cdef cppclass _ChamberParmType "ChamberParmType":
        _ChamberParmType()

    cdef cppclass _HB_ParmType "HB_ParmType":
        _HB_ParmType()

    cdef cppclass _LES_AtomType "LES_AtomType":
        _LES_AtomType()

    cdef cppclass _CmapGridType "CmapGridType":
        _CmapGridType()

    cdef cppclass _CmapType "CmapType":
        _CmapType()

# Cython wrapper classes - only include those that are actually used
cdef class AngleType:
    cdef _AngleType* thisptr

cdef class BondType:
    cdef _BondType* thisptr


cdef class DihedralType:
    cdef _DihedralType* thisptr