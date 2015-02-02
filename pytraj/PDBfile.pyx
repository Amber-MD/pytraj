# distutils: language = c++


cdef class PDBfile:
    def __cinit__(self):
        self.thisptr = new _PDBfile()

    def __dealloc__(self):
        del self.thisptr

    #def bint ID_PDB(self,CpptrajFile):

    #def PDB_RECTYPE NextRecord(self):

    #def Atom pdb_Atom(self,bint):

    #def Atom pdb_Atom(self):

    #def void pdb_XYZ(self,double *):

    #def void pdb_Box(self,double *):

    #def NameType pdb_ResName(self):

    #def int pdb_ResNum(self):

    #def PDB_RECTYPE RecType(self):

    #def void WriteTER(self,int, NameType, char, int):

    #def void WriteHET(self,int, double, double, double):

    #def void WriteATOM(self,int, double, double, double, char *, double):

    #def void WriteATOM(self, char *, int, double, double, double, char *, double):

    #def void WriteCoord(self,PDB_RECTYPE, int, NameType, NameType, char, int, double, double, double):

    #def void WriteCoord(self,PDB_RECTYPE, int, NameType, NameType, char, int, double, double, double, float, float, char *, int, bint):

    #def void WriteANISOU(self,int, NameType, NameType, char, int, int, int, int, int, int, int, char *, int):

    #def void WriteTITLE(self, string):

    #def void WriteCRYST1(self, double *):

    #def void WriteMODEL(self,int):

    #def void WriteENDMDL(self):

    #def void WriteEND(self):

