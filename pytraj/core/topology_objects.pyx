# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray
from pytraj.c_dict import get_key, AtomicElementDict
from pytraj.externals.six import string_types
from pytraj.core.elements import Element

cdef class Atom:
    '''Atom

    Examples
    --------
    >>> from pytraj.core import Atom
    >>> # name, type, charge, mass
    >>> H = Atom(name='H', type='H', charge='0.0', mass='1.0', resid=0)
    '''

    def __cinit__(self, **kwd):
        cdef NameType aname, atype

        # atom name, atom type, charge, mass
        if kwd:
            aname = NameType(kwd.get('name', ''))
            atype = NameType(kwd.get('type', ''))
            charge = kwd.get('charge', 0.)
            mass = kwd.get('mass', 0.0)
            self.thisptr = new _Atom(aname.thisptr[0], charge, mass, atype.thisptr[0])
            self.resid = kwd.get('resid', 0)
            self._index = 0
        else:
            self.thisptr = new _Atom()
            self._index = 0

        self.own_memory = True

    def __dealloc__(self):
        if self.own_memory:
            del self.thisptr

    def copy(self):
        cdef Atom atom = Atom()
        #del atom.thisptr
        atom.thisptr = new _Atom(self.thisptr[0])
        return atom

    property index:
        def __get__(self):
            return self._index

        def __set__(self, int idx):
            self._index = idx

    def bonded_indices(self):
        """get bond indices that `self` bonds to
        """
        cdef pyarray arr0 = pyarray('i', [])
        cdef bond_iterator it

        it = self.thisptr.bondbegin()
        while it != self.thisptr.bondend():
            arr0.append(deref(it))
            incr(it)
        return arr0

    property resid:
        def __set__(self, int residIn):
            self.thisptr.SetResNum(residIn)

        def __get__(self):
            return self.thisptr.ResNum()

    def set_mol(self, int mol_num):
        self.thisptr.SetMol(mol_num)

    property charge:
        def __set__(self, double qin):
            self.thisptr.SetCharge(qin)

        def __get__(self):
            return self.thisptr.Charge()

    property gb_radius:
        def __set__(self, double r):
            self.thisptr.SetGBradius(r)

        def __get__(self):
            return self.thisptr.GBRadius()

    def __str__(self):
        if self.atomic_number > 0:
            name = self.name.strip()
            if name != '':
                txt = 'Atom(name={}, type={}, atomic_number={}, index={}, resid={})'.format(
                    name, self.type.strip(), self.atomic_number, self.index, self.resid)
            else:
                txt = 'Atom()'
        else:
            txt = 'Atom()'
        return txt

    def __repr__(self):
        return str(self)

    @property
    def element(self):
        # I really miss python 3.5 unpacking
        # return name, *Element[self.atomic_number]
        name = get_key(self.thisptr.Element(), AtomicElementDict)
        return name.lower()

    @property
    def atomic_number(self):
        return self.thisptr.AtomicNumber()

    @property
    def name(self):
        name = self.thisptr.c_str().decode('UTF-8')
        return name.rstrip()

    @property
    def type(self):
        return self.thisptr.Type().Truncated().decode()

    @property
    def molnum(self):
        return self.thisptr.MolNum()

    @property
    def chain(self):
        return self.thisptr.MolNum()

    @property
    def n_bonds(self):
        return self.thisptr.Nbonds()

    @property
    def mass(self):
        return self.thisptr.Mass()

    def is_bonded_to(self, int idx):
        return self.thisptr.IsBondedTo(idx)

    @classmethod
    def _get_bond_length(cls, id1, id2):
        """get_bond_length(id1, id2)
        Return : bond length of two atomic elements (Angstrom)

        Parameters:
        ---------
        id1 : str, AtomicElement 1
        id2 : str, AtomicElement 2
        """
        id1 = id1.upper()
        id2 = id2.upper()
        return _Atom.GetBondLength(AtomicElementDict[id1],
                                   AtomicElementDict[id2])


cdef class Residue:
    '''Residue

    Examples
    --------
    >>> Residue('ALA', resid=0, icode=0, chainID=0)
    '''

    def __cinit__(self, name='', int resid=0, icode=0, chainID=0):
        cdef NameType resname = NameType(name)
        cdef char icode_ = <int> icode
        cdef char chainID_ = <int> chainID
        self.thisptr = new _Residue(resname.thisptr[0], <int> resid,
                                    icode_, chainID_)

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        if self.n_atoms > 0:
            name = self.name.split()[0]
            txt = "<%s%s, %s atoms>" % (name,
                                        self.original_resid-1,
                                        self.n_atoms)
        else:
            txt = '<Emtpy Residue>'
        return txt

    def __repr__(self):
        return self.__str__()

    def set_last_atom(self, int i):
        self.thisptr.SetLastAtom(i)

    @property
    def first_atom_index(self):
        return self.thisptr.FirstAtom()

    @property
    def last_atom_index(self):
        return self.thisptr.LastAtom()

    @property
    def first(self):
        """return first atom index (alias of first_atom_index. (experiment))
        """
        return self.thisptr.FirstAtom()

    @property
    def last(self):
        """return last atom index (alias of last_atom_index. (experiment))
        """
        return self.thisptr.LastAtom()

    property original_resid:
        def __get__(self):
            return self.thisptr.OriginalResNum()

        def __set__(self, int i):
            self.thisptr.SetOriginalNum(i)

    @property
    def index(self):
        """shortcut of original_resid"""
        return self.thisptr.OriginalResNum() - 1

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
        name = self.thisptr.c_str().decode('UTF-8')
        return name.rstrip()


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

    def set_first(self, int begin):
        self.thisptr.SetFirst(begin)

    def set_last(self, int last):
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
