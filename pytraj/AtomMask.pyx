# distutils: language = c++
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarrary

# FIXME : property does not work properly


cdef class AtomMask(object):
    # TODO : rename methods, add doc
    def __cinit__(self, *args):
        cdef int begin_atom, end_atom, atom_num
        cdef string maskstring
        cdef AtomMask rhs_atm
        if not args:
            self.thisptr = new _AtomMask()
        else:
            if len(args) == 1:
                if isinstance(args[0], int):
                    atom_num = args[0]
                    self.thisptr = new _AtomMask(atom_num)
                elif isinstance(args[0], AtomMask):
                    rhs_atm = args[0]
                    self.thisptr = new _AtomMask(rhs_atm.thisptr[0])
                else:
                    # string_types
                    maskstring = args[0].encode("UTF-8")
                    self.thisptr = new _AtomMask(maskstring)
            elif len(args) == 2:
                begin_atom, end_atom = args
                self.thisptr = new _AtomMask(begin_atom, end_atom)
            else:
                # TODO: better Error
                raise NotImplementedError()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def selected_indices(self):
        return pyarrary('i', self.thisptr.Selected())

    def __iter__(self):
        cdef cppvector[int].const_iterator it
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            yield deref(it)
            incr(it)

    def back(self):
        return self.thisptr.back()

    @property
    def n_selected(self):
        return self.thisptr.Nselected()

    #def  int operator[](self,int idx):
    def __getitem__(self, int idx):
        return self.thisptr.index_opr(idx)

    def mask_string(self):
        return self.thisptr.MaskString()

    def mask_expression(self):
        return self.thisptr.MaskExpression()

    def mask_string_set(self):
        return self.thisptr.MaskStringSet()

    def none(self):
        return self.thisptr.None()

    def is_char_mask(self):
        return self.thisptr.IsCharMask()

    def reset_mask(self):
        self.thisptr.ResetMask()

    def clear_selected(self):
        self.thisptr.ClearSelected()

    def invert_mask(self):
        self.thisptr.InvertMask()

    def num_atoms_in_common(self, AtomMask other_mask):
        return self.thisptr.NumAtomsInCommon(other_mask.thisptr[0])

    def add_selected_atom(self,int i):
        self.thisptr.AddSelectedAtom(i)

    def add_atom(self,int atom_num):
        self.thisptr.AddAtom(atom_num)

    def add_atoms(self, vector[int] v):
        self.thisptr.AddAtoms(v)

    def add_atom_range(self, int begin, int end):
        self.thisptr.AddAtomRange(begin, end)

    def add_mask_at_position(self, AtomMask atm, int pos):
        self.thisptr.AddMaskAtPosition(atm.thisptr[0], pos)

    def print_mask_atoms(self, char* mask):
        self.thisptr.PrintMaskAtoms(mask)

    def set_mask_string(self, mask):
        mask = mask.encode()
        return self.thisptr.SetMaskString(mask)

    def setup_int_mask(self, char *charmask, int natom, int debug=0):
        self.thisptr.SetupIntMask(charmask, natom, debug)

    def setup_char_mask(self, char* charmask, int natom, int debug=0):
        self.thisptr.SetupCharMask(charmask, natom, debug)

    def atoms_in_char_mask(self, *args):
        cdef atomid, begin, end
        if len(args) == 1:
            atomid = args[0]
            return self.thisptr.AtomInCharMask(atomid)
        elif len(args) == 2:
            begin, end = args
            return self.thisptr.AtomsInCharMask(begin, end)

    def set_natom(self,int a):
        self.thisptr.SetNatom(a)

    def convert_to_char_mask(self):
        return self.thisptr.ConvertToCharMask()

    def convert_to_int_mask(self):
        return self.thisptr.ConvertToIntMask()

    def mask_info(self):
        self.thisptr.MaskInfo()

    def brief_mask_info(self):
        self.thisptr.BriefMaskInfo()

    # Not yet support
    #def  token_iterator begintoken(self):
    #def  token_iterator endtoken(self):

