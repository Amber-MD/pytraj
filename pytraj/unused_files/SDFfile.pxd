# distutils: language = c++
from CpptrajFile cimport *
from Atom cimport *


cdef extern from "SDFfile.h": 
    cdef cppclass _SDFfile "SDFfile":
        _SDFfile() 
        bint ID_SDF(_CpptrajFile &)
        bint ReadHeader() 
        int SDF_XYZ(double *)
        _Atom SDF_Atom() 
        int SDF_Bond(int &, int &)
        int SDF_Natoms() const 
        int SDF_Nbonds() const 
        const string& SDF_Title() const 
