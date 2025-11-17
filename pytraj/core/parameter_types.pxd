# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string


cdef extern from "ParameterTypes.h":
    # Classes used in topology.pyx
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

    # Array classes with iterator support
    cdef cppclass _BondArray "BondArray":
        _BondArray()
        ctypedef _BondType* iterator
        ctypedef const _BondType* const_iterator
        iterator begin()
        iterator end()
        const_iterator begin() const
        const_iterator end() const
        void insert(iterator pos, const_iterator first, const_iterator last)

    cdef cppclass _AngleArray "AngleArray":
        _AngleArray()
        ctypedef _AngleType* iterator
        ctypedef const _AngleType* const_iterator
        iterator begin()
        iterator end()
        const_iterator begin() const
        const_iterator end() const
        void insert(iterator pos, const_iterator first, const_iterator last)

    cdef cppclass _DihedralArray "DihedralArray":
        _DihedralArray()
        ctypedef _DihedralType* iterator
        ctypedef const _DihedralType* const_iterator
        iterator begin()
        iterator end()
        const_iterator begin() const
        const_iterator end() const

    # Minimal classes for completeness but not directly used in topology.pyx
    cdef cppclass _BondParmType "BondParmType":
        _BondParmType()

    cdef cppclass _AngleParmType "AngleParmType":
        _AngleParmType()

    cdef cppclass _NonbondType "NonbondType":
        _NonbondType()

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