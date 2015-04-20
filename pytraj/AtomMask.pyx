# distutils: language = c++
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarrary
from cpython cimport array
from pytraj.decorators import deprecated
from pytraj._utils import set_world_silent

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

    @property
    def indices(self):
        return pyarrary('i', self.thisptr.Selected())

    @property
    def _indices_view(self):
        cdef vector[int] v = self.thisptr.Selected()
        cdef pyarrary a_empty = pyarrary('i', [])
        cdef int size = v.size()
        cdef pyarrary arr0 = array.clone(a_empty, size, zero=True) 
        cdef int[:] myview = arr0
        cdef int i 

        for i in range(size):
            myview[i] = v[i]

        return myview

    def __iter__(self):
        cdef cppvector[int].const_iterator it
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            yield deref(it)
            incr(it)

    def back(self):
        return self.thisptr.back()

    @property
    @deprecated
    def n_selected(self):
        return self.thisptr.Nselected()

    property n_atoms:
        def __get__(self):
            """the number of selected atoms based on mask"""
            return self.thisptr.Nselected()
        def __set__(self, int value):
            self.thisptr.SetNatom(value)

    def __getitem__(self, int idx):
        return self.thisptr.index_opr(idx)

    property mask_string:
        def __get__(self):
            return self.thisptr.MaskString()
        def __set__(self, value):
            self.thisptr.SetMaskString(value.encode())

    def mask_expression(self):
        cdef string t
        t = self.thisptr.MaskExpression()
        return t.decode()

    def mask_string_was_set(self):
        return self.thisptr.MaskStringSet()

    def is_empty(self):
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

    def add_selected_indices(self, arr0):
        """add atom index without sorting

        See Also
        --------
        add_atom
        """
        cdef int[:] int_view
        cdef int i

        try:
            # try casting to memview
            int_view = arr0
            for i in range(int_view.shape[0]):
                self.thisptr.AddSelectedAtom(int_view[i])
        except TypeError:
            # slower way if array does not have buffer interface
            for i in arr0:
                self.thisptr.AddSelectedAtom(i)

    def add_atom(self,int atom_num):
        """add atom index and sort"""
        self.thisptr.AddAtom(atom_num)

    def add_atoms(self, vector[int] v):
        self.thisptr.AddAtoms(v)

    def add_atom_range(self, int begin, int end):
        self.thisptr.AddAtomRange(begin, end)

    def add_mask_at_position(self, AtomMask atm, int pos):
        self.thisptr.AddMaskAtPosition(atm.thisptr[0], pos)

    def print_mask_atoms(self, mask):
        set_world_silent(False)
        mask = mask.encode()
        self.thisptr.PrintMaskAtoms(mask)
        set_world_silent(True)

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

    def convert_to_char_mask(self):
        return self.thisptr.ConvertToCharMask()

    def convert_to_int_mask(self):
        return self.thisptr.ConvertToIntMask()

    def mask_info(self):
        set_world_silent(False)
        self.thisptr.MaskInfo()
        set_world_silent(True)

    def brief_mask_info(self):
        set_world_silent(False)
        self.thisptr.BriefMaskInfo()
        set_world_silent(True)
