# distutils: language = c++
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarrary
from cpython cimport array
from pytraj.decorators import deprecated
from pytraj._set_silent import set_world_silent
from pytraj.externals.six import string_types
from pytraj.utils import is_array
from pytraj._utils import _int_array1d_like_to_memview
from pytraj.compat import range

__all__ = ['AtomMask']

cdef class AtomMask(object):
    def __cinit__(self, *args):
        cdef int begin_atom, end_atom, atom_num
        cdef string maskstring
        cdef AtomMask rhs_atm

        if not args:
            self.thisptr = new _AtomMask()
        elif is_array(args[0]) or isinstance(args[0], (list, tuple, range)):
            self.thisptr = new _AtomMask()
            self.add_selected_indices(args[0])
        else:
            if len(args) == 1:
                if isinstance(args[0], int):
                    atom_num = args[0]
                    self.thisptr = new _AtomMask(atom_num)
                elif isinstance(args[0], AtomMask):
                    rhs_atm = args[0]
                    self.thisptr = new _AtomMask(rhs_atm.thisptr[0])
                elif isinstance(args[0], string_types):
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

    @property
    @deprecated
    def n_selected(self):
        return self.thisptr.Nselected()

    property n_atoms:
        def __get__(self):
            """the number of selected atoms based on mask"""
            return self.thisptr.Nselected()

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

    def is_empty(self):
        return self.thisptr.None()

    def invert_mask(self):
        self.thisptr.InvertMask()

    def add_selected_indices(self, arr0):
        """add atom index without sorting

        See Also
        --------
        add_atom
        """
        cdef int[:] int_view
        cdef int i

        # try casting to memview
        if not is_array(arr0):
            int_view = _int_array1d_like_to_memview(arr0)
        else:
            try:
                int_view = arr0
            except:
                # numpy compat
                int_view = arr0.astype('i4')

        for i in range(int_view.shape[0]):
            self.thisptr.AddSelectedAtom(int_view[i])

    def add_atom(self,int atom_num):
        """add atom index and sort"""
        self.thisptr.AddAtom(atom_num)

    def add_atoms(self, vector[int] v):
        self.thisptr.AddAtoms(v)

    def add_atom_range(self, int begin, int end):
        self.thisptr.AddAtomRange(begin, end)

    def add_mask_at_position(self, AtomMask atm, int pos):
        self.thisptr.AddMaskAtPosition(atm.thisptr[0], pos)
