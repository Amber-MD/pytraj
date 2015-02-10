# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from libcpp.string cimport string
from cpython.array cimport array as pyarray
#from pytraj.TopologyList cimport TopologyList

from pytraj.decorators import name_will_be_changed
from pytraj.utils.check_and_assert import _import_numpy
from pytraj._ParmFile import TMPParmFile
try:
    set
except NameError:
    from sets import Set as set
from pytraj.externals.six import PY3, PY2, string_types, binary_type

cdef class Topology:
    def __cinit__(self, *args):
        """
        args = filename or Topology instance
        """
        #cdef TopologyList toplist = TopologyList()
        cdef Topology tp
        cdef string filename
        self.thisptr = new _Topology()

        # I dont make default py_free_mem (True) in __cinit__ since
        # when passing something like top = Topology(filename), Python/Cython
        # would think "filename" is value of py_free_mem
        self.py_free_mem = True

        if not args:
            #print "there is no args" # for debug
            # make empty Topology instance
            pass
        else:
            if len(args) == 1:
                if isinstance(args[0], string_types):
                    filename = args[0].encode()
                    pf = TMPParmFile()
                    tp = Topology()
                    pf.readparm(filename, tp)
                    self.thisptr[0] = tp.thisptr[0]
                elif isinstance(args[0], Topology):
                    tp = args[0]
                    self.thisptr[0] =  tp.thisptr[0]
            else:
                raise ValueError()

    def __dealloc__(self):
        if self.py_free_mem and self.thisptr:
            del self.thisptr

    def __str__(self):
        tmp = "%s instance with %s atoms. ID = %s" % (
                self.__class__.__name__,
                self.n_atoms,
                hex(id(self))
                )
        return tmp

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.n_atoms

    def load(self, string filename):
        """loading Topology from filename

        filename : {str}

        if Topology instance is not empty, it will be still replaced by new one

        # seriously why do we need this method?
        >>> top = Topology("./data/Tc5b.top")
        >>> # replace old with new topology
        >>> top.load("./data/HP36.top")
        >>> # why not using "top = Topology("./data/HP36.top")"?
        """
        self = Topology(filename)

    def copy(self, *args):
        """return a copy of 'self' or copy from 'other' to 'self'
        TODO : add more doc
        """
        cdef Topology tmp
        cdef Topology other

        if not args:
            # make a copy of 'self'
            tmp = Topology()
            tmp.thisptr[0] = self.thisptr[0]
            return tmp
        elif isinstance(args[0], Topology):
            # copy other Topology instance to "self"
            # no error? really?
            other = args[0]
            self.thisptr[0] = other.thisptr[0]

    def __getitem__(self, idx):
        """
        return Atom instance

        TODO : return either atoms or residues 

                Example : 
                    self['atom'], self['residue'] 
                        return list of Atom instances or Residue instances

                    self[0] 
                        return 0-th Atom instance
        END TODO
        """

        cdef Atom atom 
        cdef int i

        if isinstance(idx, (int, long)):
            # need to explicitly cast to int
            i = <int> idx
            atom = Atom()
            atom.thisptr[0] = self.thisptr.index_opr(i)
            return atom
        elif isinstance(idx, string_types):
            alist = []
            # return atom object iterator with given mask
            # self(idx) return AtomMask object
            for i in self(idx).selected_indices():
                alist.append(self[i])
            return alist
        else:
            raise ValueError("must be integer or string")

    def __call__(self, mask, *args, **kwd):
        """intended to use with Frame indexing
        Return : AtomMask object
        >>> frame[top("@CA")]
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        return atm

    def atomiter(self):
        cdef Atom atom
        cdef atom_iterator it

        it = self.thisptr.begin()
        while it != self.thisptr.end():
            atom = Atom()
            atom.thisptr[0] = deref(it)
            yield atom
            incr(it)

    def residueiter(self):
        cdef Residue res
        cdef res_iterator it
        it = self.thisptr.ResStart()

        while it != self.thisptr.ResEnd():
            res = Residue()
            res.thisptr[0] = deref(it)
            yield res
            incr(it)
        
    @property
    def moliter(self):
        cdef Molecule mol
        cdef mol_iterator it
        it = self.thisptr.MolStart()

        while it != self.thisptr.MolEnd():
            mol = Molecule()
            mol.thisptr[0] = deref(it)
            yield mol
            incr(it)

    def set_parm_name(self, string title, FileName filename):
        self.thisptr.SetParmName(title, filename.thisptr[0])

    def set_reference_coords(self, Frame frameIn):
        self.thisptr.SetReferenceCoords(frameIn.thisptr[0])

    def file_path(self):
        return self.thisptr.c_str()

    def trunc_res_atom_name(self, int atom):
        return self.thisptr.TruncResAtomName(atom)

    def atom_mask_name(self, int atom):
        return self.thisptr.AtomMaskName(atom)

    def trunc_resname_num(self, int res):
        return self.thisptr.TruncResNameNum(res)

    def find_atom_in_residue(self, int res, NameType atname):
        return self.thisptr.FindAtomInResidue(res, atname.thisptr[0])
    
    def find_residue_max_natom(self):
        return self.thisptr.FindResidueMaxNatom()
    
    #def SoluteAtoms(self):
    #    return self.thisptr.SoluteAtoms()

    @property
    def atomlist(self):
        """return list of atoms
        """
        cdef Atom atom 
        cdef vector[_Atom].iterator it 
        cdef vector[_Atom] v
        cdef list atlist = []

        v = self.thisptr.Atoms()
        it = v.begin()
        while it != v.end():
            atom = Atom()
            atom.thisptr[0] = deref(it)
            atlist.append(atom)
            incr(it)
        return atlist

    @property
    def residuelist(self):
        reslist = []
        for res in self.residueiter():
            reslist.append(res)
        return reslist

    @property
    def moleculelist(self):
        mlist = []
        for mol in self.moliter:
            mlist.append(mol)
        return mlist

    def summary(self):
        self.thisptr.Summary()

    def brief(self, char* heading="*"):
        self.thisptr.Brief(heading)

    def atom_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintAtomInfo(maskString)

    def bond_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintBondInfo(maskString)
    
    def angle_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintAngleInfo(maskString)

    def dihedral_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintDihedralInfo(maskString)

    def molecule_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintMoleculeInfo(maskString)

    def residue_info(self, maskString="*"):
        maskString = maskString.encode()
        self.thisptr.PrintResidueInfo(maskString)

    def charge_mass_info(self, maskString, int idtype):
        maskString = maskString.encode()
        self.thisptr.PrintChargeMassInfo(maskString, idtype)

    # BROKEN
    #def has_vel(self):
    #    return self.thisptr.HasVelInfo()
    
    def add_atom(self, Atom atom=Atom(), 
                 int resid=0, 
                 resname="", xyz=None):
        """add_atom(Atom atomIn, int o_resnum, NameType resname, double[:] XYZin)"""
        cdef double[:] XYZin
        cdef NameType restype = NameType(resname)

        if xyz is None:
            XYZin = pyarray('d', [0., 0., 0.])
        else:
            XYZin = pyarray('d', xyz)

        self.thisptr.AddTopAtom(atom.thisptr[0], resid, restype.thisptr[0], &XYZin[0])

    def start_new_mol(self):
        self.thisptr.StartNewMol()

    def common_setup(self, bint bondsearch):
        return self.thisptr.CommonSetup(bondsearch)

    def set_offset(self, double x):
        self.thisptr.SetOffset(x)

    def set_ipol(self, int id):
        self.thisptr.SetIpol(id)

    def orig_filename(self):
        cdef FileName filename = FileName()
        filename.thisptr[0] = self.thisptr.OriginalFilename()
        return filename

    property parm_index:
        def __get__(self):
            return self.thisptr.Pindex()

    property p_index:
        # shortcut of parm_index
        def __get__(self):
            return self.thisptr.Pindex()

    property n_atoms:
        def __get__(self):
            return self.thisptr.Natom()

    property n_res:
        # shortcur 
        def __get__(self):
            return self.n_residues

    property n_residues:
        def __get__(self):
            return self.thisptr.Nres()

    property n_mols:
        def __get__(self):
            return self.thisptr.Nmol()

    property n_solvents:
        def __get__(self):
            return self.thisptr.Nsolvent()

    property n_frames:
        def __get__(self):
            return self.thisptr.Nframes()

    # BROKEN
    #property n_repdims:
    #    def __get__(self):
    #        return self.thisptr.NrepDims()

    property parm_name:
        def __get__(self):
            return self.thisptr.ParmName()
        def __set__(self, name):
            # TODO : check
            self.thisptr.SetParmName(name, FileName().thisptr[0])

    property GB_radiiset:
        def __get__(self):
            return self.thisptr.GBradiiSet()

    #def int SetAmberExtra(self, vector[double], vector[NameType], vector[int], vector[int]):

    def set_integer_mask(self, AtomMask atm, Frame frame=Frame()):
        if frame.is_empty():
            return self.thisptr.SetupIntegerMask(atm.thisptr[0])
        else:
            return self.thisptr.SetupIntegerMask(atm.thisptr[0], frame.thisptr[0])

    def set_char_mask(self, AtomMask atm, Frame frame=Frame()):
        if frame.is_empty():
            return self.thisptr.SetupCharMask(atm.thisptr[0])
        else:
            return self.thisptr.SetupCharMask(atm.thisptr[0], frame.thisptr[0])

    def scale_dihedral_k(self, double value):
        self.thisptr.ScaleDihedralK(value)

    def set_box(self, Box boxin):
        self.thisptr.SetParmBox(boxin.thisptr[0])

    def partial_modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.partialModifyStateByMask(m.thisptr[0]))
        return top

    def modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.modifyStateByMask(m.thisptr[0]))
        return top

    def modify_by_map(self, vector[int] m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.ModifyByMap(m))
        return top

    def strip_atoms(Topology self, mask, copy=False):
        # TODO : shorter way?
        """strip atoms with given mask"""
        cdef AtomMask atm = AtomMask()
        cdef Topology tmptop = Topology()
        mask = mask.encode()

        atm.thisptr.SetMaskString(mask)
        atm.thisptr.InvertMask()
        self.thisptr.SetupIntegerMask(atm.thisptr[0])
        tmptop.thisptr = self.thisptr.modifyStateByMask(atm.thisptr[0])
        if copy:
            return tmptop
        else:
            self.thisptr[0] = tmptop.thisptr[0]

    def tag(self):
        # what does this do?
        return self.thisptr.Tag()

    def is_empty(self):
        s = self.file_path()
        return s == ""

    def get_atom_indices(self, mask, *args, **kwd):
        """return atom indices with given mask
        To be the same as cpptraj/Ambertools: we mask indexing starts from 1
        but the return list/array use 0
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        has_numpy, np = _import_numpy()
        if has_numpy:
            # ndarray
            return np.asarray(atm.selected_indices())
        else:
            # list
            return atm.selected_indices()

    @name_will_be_changed("")
    def get_unique_resname(self):
        s = set()
        for res in self.residueiter():
            s.add(res.name)
        return s

    @name_will_be_changed("")
    def get_unique_atomname(self):
        s = set()
        for atom in self.atomiter():
            s.add(atom.name)
        return s

    def get_atomname_set(self):
        return self.get_unique_atomname()

    def get_resname_set(self):
        return self.get_unique_resname()

    def parm_coordinnfo(self):
        cdef CoordinateInfo coordinfo = CoordinateInfo()
        coordinfo.thisptr[0] = self.thisptr.ParmCoordInfo()
        return coordinfo

    def join(self, top):
        cdef Topology _top
        if isinstance(top, Topology):
            _top = top
            if _top == self:
                raise ValueError("can not join yourself, use copy() method")
        elif isinstance(top, string_types):
            _top = Topology(top)
        else:
            raise ValueError("support only Topology object or top filename")

        self.thisptr.AppendTop(_top.thisptr[0])
