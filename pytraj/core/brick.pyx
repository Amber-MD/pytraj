# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray
from pytraj.cpptraj_dict import get_key, AtomicElementDict
from pytraj.externals.six import string_types

cdef class Atom:
#    def __cinit__(self, aname="", atype=None, index=0, 
#                        chainID=' ', resnum=0, mol=0,
#                        charge=0.0, polar=0.0, mass=1.0, 
#                        element=None,
#                        *args, **kwd):
    def __cinit__(self, *args, **kwd):
        # TODO : add more constructors
        cdef NameType aname, atype
        cdef Atom other
        cdef double charge, mass

        if not args and not kwd:
            self.thisptr = new _Atom()
        else:
            if len(args) == 1 and isinstance(args[0], Atom):
                other = <Atom>  args[0]
                self.thisptr = new _Atom(other.thisptr[0])
            elif len(args) == 2:
                if isinstance(args[0], string_types) and isinstance(args[1], string_types):
                    aname = NameType(args[0])
                    atype = NameType (args[1])
                    self.thisptr = new _Atom(aname.thisptr[0], atype.thisptr[0], 1.0) 
            elif len(args) == 4:
                # atom name, atom type, charge, mass
                aname = NameType(args[0])
                atype = NameType (args[1])
                charge = args[2]
                mass = args[3]
                self.thisptr = new _Atom(aname.thisptr[0], charge, mass, atype.thisptr[0])
            else:
                raise NotImplementedError("not yet supported")

    def __dealloc__(self):
        del self.thisptr

    def copy(self):
        cdef Atom atom = Atom()
        #del atom.thisptr
        atom.thisptr = new _Atom(self.thisptr[0])
        return atom

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

    property resnum:
        def __set__(self, int resnumIn):
            self.thisptr.SetResNum(resnumIn)
        def __get__(self):
            return self.thisptr.ResNum()

    property mol:
        def __set__(self,int molIn):
            self.thisptr.SetMol(molIn)

    property charge:
        def __set__(self,double qin):
            self.thisptr.SetCharge(qin)
        def __get__(self):
            return self.thisptr.Charge()
    
    def no_mol(self):
        return self.thisptr.NoMol()

    def __str__(self):
        if self.atomic_number > 0:
            name = self.name
            if name != '':
                name = name.split()[0]
                txt = "<%s-atom, resnum=%s, n_bonds=%s>" % (name, self.resnum, self.n_bonds)
            else:
                txt = '<Empty Atom>'
        else:
            txt = '<Empty Atom>'
        return txt

    def __repr__(self):
        return self.__str__()

    @property
    def element(self):
        return get_key(self.thisptr.Element(), AtomicElementDict)

    @property
    def atomic_number(self):
        return self.thisptr.AtomicNumber()

    def nametype(self):
        # TODO : do we need this method?
        cdef NameType nt = NameType()
        nt.thisptr[0] = self.thisptr.Name()
        return nt

    @property
    def name(self):
        name = self.thisptr.c_str().decode('UTF-8')
        return name.rstrip()

    @property
    def type(self):
        cdef NameType nt = NameType()
        nt.thisptr[0] = self.thisptr.Type()
        return nt

    @property
    def typeindex(self):
        return self.thisptr.TypeIndex()

    @property
    def molnum(self):
        return self.thisptr.MolNum()

    @property
    def n_bonds(self):
        return self.thisptr.Nbonds()

    @property
    def mass(self):
        return self.thisptr.Mass()

    def is_bonded_to(self, int idx):
        # TODO : add doc
        return self.thisptr.IsBondedTo(idx)

    @classmethod
    def get_bond_length(cls, id1, id2):
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
# distutils: language = c++
from pytraj.externals.six import string_types


cdef class Residue:
    def __cinit__(self):
        self.thisptr = new _Residue()

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        if self.n_atoms > 0:
            name = self.name.split()[0]
            txt = "<%s%s, %s atoms>" % (name,
                                       self.original_resnum-1,
                                       self.n_atoms)
        else:
            txt = '<Emtpy Residue>'
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

