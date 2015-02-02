# distutils: language = c++

cdef class Molecule:
    def __cinit__(self, *args):
        cdef int beginidx, endidx
        if not args:
            self.thisptr = new _Molecule()
        elif len(args) == 2:
            beginidx, endidx = args
            self.thisptr = new _Molecule(beginidx, endidx)
        else:
            raise ValueError("Must have None or 2 integers")

    def __dealloc__(self):
        del self.thisptr

    def set_first(self,int begin):
        self.thisptr.SetFirst(begin)

    def set_last(self,int last):
        self.thisptr.SetLast(last)

    def set_solvent(self):
        self.thisptr.SetSolvent()

    def set_no_solvent(self):
        self.thisptr.SetNoSolvent()

    @property
    def begin_atom(self):
        return self.thisptr.BeginAtom()

    @property
    def end_atom(self):
        return self.thisptr.EndAtom()

    def is_solvent(self):
        return self.thisptr.IsSolvent()

    @property
    def n_atoms(self):
        return self.thisptr.NumAtoms()

