# distutils: language = c++
from pytraj.externals.six import string_types


cdef class Residue:
    def __cinit__(self, *args):
        cdef int orig_num
        cdef NameType nametype 
        cdef int firstAtomIn

        if not args:
            self.thisptr = new _Residue()
        else:
            if len(args) == 3:
                orig_num, _resname, firstAtomIn = args
                if isinstance(_resname, string_types):
                    #  try to cast to NameType
                    nametype = NameType(_resname)
                else:
                    nametype = _resname
                self.thisptr = new _Residue(orig_num, nametype.thisptr[0], firstAtomIn)

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        txt = "%s-%s" % (self.thisptr.c_str().decode('UTF-8'), self.original_resnum-1)
        return txt

    def __repr__(self):
        return self.__str__()

    def set_last_atom(self,int i):
        self.thisptr.SetLastAtom(i)

    def set_original_num(self, int i):
        self.thisptr.SetOriginalNum(i)

    @property
    def first_atom_idx(self):
        return self.thisptr.FirstAtom()

    @property
    def last_atom_idx(self):
        return self.thisptr.LastAtom()

    @property
    def original_resnum(self):
        return self.thisptr.OriginalResNum()

    @property
    def index(self):
        """shortcut of original_resnum"""
        return self.original_resnum()

    def ntype(self):
        cdef NameType nt = NameType()
        nt.thisptr[0] = self.thisptr.Name()
        return nt

    @property
    def n_atoms(self):
        return self.thisptr.NumAtoms()

    def is_solvent(self):
        return self.thisptr.NameIsSolvent()

    @property
    def name(self):
        return self.thisptr.c_str().decode('UTF-8')
