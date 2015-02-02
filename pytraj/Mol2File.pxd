# distutils: language = c++
from CpptrajFile cimport *
from Atom cimport *


cdef extern from "Mol2File.h": 
    # Mol2File.h
    ctypedef enum TRIPOSTAG "Mol2File::TRIPOSTAG":
        MOLECULE "Mol2File::MOLECULE"
        ATOM "Mol2File::ATOM"
        BOND "Mol2File::BOND"
        SUBSTRUCT "Mol2File::SUBSTRUCT"
    cdef cppclass _Mol2File "Mol2File":
        _Mol2File() 
        bint IsMol2Keyword(const char *)
        bint ID_Mol2(_CpptrajFile&)
        int ScanTo(TRIPOSTAG)
        void WriteHeader(TRIPOSTAG)
        bint Read_Molecule() 
        bint Write_Molecule(bint, int)
        int Next_Molecule() 
        int Mol2Bond(int&, int&)
        int Mol2XYZ(double *)
        _Atom Mol2_Atom() 
        _NameType Mol2_Residue(int&)
        void SetMol2Natoms(int nIn)
        void SetMol2Nbonds(int nIn)
        void SetMol2Title(const string& tIn)
        int Mol2Natoms() const 
        int Mol2Nbonds() const 
        const string& Mol2Title() const 
