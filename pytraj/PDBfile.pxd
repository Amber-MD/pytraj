# distutils: language = c++
from libcpp.string cimport string
from pytraj.CpptrajFile cimport *
from pytraj.Atom cimport *


cdef extern from "PDBfile.h": 
    # PDBfile.h
    ctypedef enum PDB_RECTYPE "PDBfile::PDB_RECTYPE":
        pass
        #ATOM "PDBfile::ATOM"
        #HETATM "PDBfile::HETATM"
        #CRYST1 "PDBfile::CRYST1"
        #TER "PDBfile::TER"
        #END "PDBfile::END"
        #ANISOU "PDBfile::ANISOU"
        #END_OF_FILE "PDBfile::END_OF_FILE"
        #UNKNOWN "PDBfile::UNKNOWN"
    cdef cppclass _PDBfile "PDBfile":
        _PDBfile() 
        bint ID_PDB(_CpptrajFile&)
        PDB_RECTYPE NextRecord() 
        Atom pdbAtom(bint)
        Atom pdbAtom() 
        void pdbXYZ(double *)
        void pdbBox(double *) const 
        _NameType pdbResName() 
        int pdb_ResNum() 
        PDB_RECTYPE RecType() const 
        void WriteTER(int, const _NameType&, char, int)
        void WriteHET(int, double, double, double)
        void WriteATOM(int, double, double, double, const char *, double)
        void WriteATOM(const char *, int, double, double, double, const char *, double)
        void WriteCoord(PDB_RECTYPE, int, const _NameType&, const _NameType&, char, int, double, double, double)
        void WriteCoord(PDB_RECTYPE, int, const _NameType&, const _NameType&, char, int, double, double, double, float, float, const char *, int, bint)
        void WriteANISOU(int, const _NameType&, const _NameType&, char, int, int, int, int, int, int, int, const char *, int)
        void WriteTITLE(const string&)
        void WriteCRYST1(const double *)
        void WriteMODEL(int)
        void WriteENDMDL() 
        void WriteEND() 


cdef class PDBfile:
    cdef _PDBfile* thisptr
