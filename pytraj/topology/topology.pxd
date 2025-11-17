# distutils: language = c++

from ..cython_extra_header.cpp_vector cimport vector as cppvector
from ..core.topology_objects cimport _Atom, Atom, _Residue, Residue, _Molecule, Molecule
from ..core.box cimport _Box, Box
from ..core.parameter_types cimport *
from ..core.c_core cimport (_FileName, FileName, _NameType, NameType)
from ..core.c_core cimport _AtomMask, AtomMask
from ..trajectory.frame cimport _Frame, Frame
from libcpp.string cimport string
from ..core.c_core cimport _FileName, FileName, _ArgList, ArgList


ctypedef cppvector[_Atom].const_iterator atom_iterator
ctypedef cppvector[_Residue].const_iterator res_iterator
ctypedef cppvector[_Molecule].const_iterator mol_iterator
ctypedef cppvector[int].const_iterator bond_iterator

cdef extern from "CoordinateInfo.h":
    cdef cppclass _CoordinateInfo "CoordinateInfo" nogil:
        _CoordinateInfo()
        _CoordinateInfo(const _Box& b, bint v, bint t, bint m)
        bint HasBox() const
        const _Box& TrajBox() const
        bint HasVel() const
        bint HasTemp() const
        bint HasTime() const
        bint HasForce() const
        bint HasReplicaDims() const
        void SetTime(bint m)
        void SetTemperature(bint t)
        void SetVelocity(bint v)
        void SetBox(const _Box& b)

cdef extern from "Topology.h":
    cdef cppclass _Topology "Topology" nogil:
        _Topology()
        int DetermineMolecules()
        void SetDistMaskRef(_Frame)
        _Atom& GetAtomView "SetAtom" (int idx)
        int Natom() const
        int Nres() const
        int Nmol() const
        int Nsolvent() const
        const _FileName& OriginalFilename() const
        atom_iterator begin()
        atom_iterator end()
        const _Atom& index_opr "operator[]"(int idx)
        const vector[_Atom]& Atoms()
        inline res_iterator ResStart()
        inline res_iterator ResEnd()
        const _Residue& Res(int idx)
        inline mol_iterator MolStart() const
        inline mol_iterator MolEnd() const
        const _BondArray& Bonds() const
        const _BondArray& BondsH() const
        void AddBond(int, int)
        const _AngleArray& Angles() const
        const _AngleArray& AnglesH() const
        const _DihedralArray& Dihedrals() const
        const _DihedralArray& DihedralsH() const
        void Summary() const
        inline const _Box& ParmBox() const
        void SetParmBox(_Box& bIn)
        int AddTopAtom(_Atom&, _Residue&)
        void AddAngle(int, int, int)
        void AddDihedral(int, int, int, int)
        bint SetupIntegerMask(_AtomMask&) const
        bint SetupIntegerMask(_AtomMask&, const _Frame&) const
        _Topology* partialModifyStateByMask(const _AtomMask& m) const
        _Topology* modifyStateByMask(const _AtomMask& m) const
        int AppendTop(const _Topology &)
        int SetSolvent(const string&)

cdef class Topology:
    cdef _Topology* thisptr
    cdef public bint _own_memory
    cdef cppvector[int] _get_atom_bond_indices(self, _Atom)

cdef extern from "ParmFile.h":
    ctypedef enum ParmFormatType "ParmFile::ParmFormatType":
        pass
        UNKNOWN_PARM "ParmFile::UNKNOWN_PARM"
    cdef cppclass _ParmFile "ParmFile" nogil:
        _ParmFile()
        int ReadTopology(_Topology&, const string&, const _ArgList&, int)
        int ReadTopology(_Topology& t, const string& n, int d)
        int WriteTopology(const _Topology&, const string&, const _ArgList&, ParmFormatType, int)
        const _FileName ParmFilename()

cdef class ParmFile:
    cdef _ParmFile* thisptr