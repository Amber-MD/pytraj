# distutils: language = c++


cdef class AtomMap:
    def __cinit__(self):
        self.thisptr = new _AtomMap()

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, int idx):
        """TODO : support fancy indexing?"""
        pass

    def __setitem__(self, int idx, value):
        pass
    #def MapAtom operator[](self,int):

    #def  MapAtom operator[](self,int i):

    @property
    def n_atoms(self):
        return self.thisptr.Natom()

    #def int Setup(self, Topology):

    #def int SetupResidue(self, Topology, int):

    #def void ResetMapping(self):
    def reset(self):
        self.thisptr.ResetMapping()

    #def bint BondIsRepeated(self,int, int):

    #def void DetermineAtomIDs(self):

    #def void MarkAtomComplete(self,int, bint):

    #def void CheckForCompleteAtoms(self):

    #def int SymmetricAtoms(self, Topology, AtomIndexArray, int):
