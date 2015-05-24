# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function
cimport cython
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from libcpp.string cimport string
from cpython.array cimport array as pyarray
from cpython cimport array as pyarray_master
from pytraj._set_silent import set_world_silent # turn on and off cpptraj's stdout
#from pytraj.TopologyList cimport TopologyList

from pytraj.decorators import name_will_be_changed
from pytraj.utils.check_and_assert import _import_numpy
from pytraj.utils.check_and_assert import is_int, is_array
from pytraj.parms._ParmFile import TMPParmFile
from pytraj.externals.six import PY3, PY2, string_types, binary_type
from pytraj.compat import set


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
        box = self.box
        if box.has_box():
            box_txt = "PBC with box type = %s" % box.type
        else:
            box_txt = "non-PBC"
         
        tmp = "<%s with %s mols, %s residues, %s atoms, %s bonds, %s>" % (
                self.__class__.__name__,
                self.n_mols,
                self.n_residues,
                self.n_atoms,
                list(self.bonds).__len__(),
                box_txt)
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
        cdef AtomMask atm
        cdef int i
        cdef list alist = []

        if is_int(idx):
            # need to explicitly cast to int
            i = <int> idx
            atom = Atom()
            atom.thisptr[0] = self.thisptr.index_opr(i)
            return atom
        elif isinstance(idx, string_types):
            # return atom object iterator with given mask
            # self(idx) return AtomMask object
            alist = [self[i] for i in self(idx)._indices_view]
        elif isinstance(idx, AtomMask):
            atm = <AtomMask> idx
            # return atom object iterator with given mask
            # self(idx) return AtomMask object
            alist = [self[i] for i in atm._indices_view]
        elif isinstance(idx, (list, tuple)) or is_array(idx):
            alist = [self[i] for i in idx]
        elif isinstance(idx, slice):
            # does not have memory efficiency with large Topology
            # (since we convert to atom list first)
            alist = self.atomlist[idx]
        if len(alist) == 1:
            return alist[0]
        else:
            return alist

    def __call__(self, mask, *args, **kwd):
        """intended to use with Frame indexing
        Return : AtomMask object
        >>> frame[top("@CA")]
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        return atm

    def __iter__(self):
        return self.atom_iter()

    @property
    def atoms(self):
        """atoms iterator, meant to be kept the same as
        other packages (ParmEd, mdtraj, ...)
        """
        return self.atom_iter()

    @property
    def residues(self):
        return self.residue_iter()

    @property
    def mols(self):
        return self.mol_iter()

    def select(self, mask):
        """return array of indices of selected atoms with `mask`

        Notes
        -----
            support openmp for distance-based atommask selction
        """
        return self(mask).indices

    def atom_iter(self):
        cdef Atom atom
        cdef atom_iterator it

        it = self.thisptr.begin()
        while it != self.thisptr.end():
            atom = Atom()
            atom.thisptr[0] = deref(it)
            yield atom
            incr(it)

    def residue_iter(self):
        cdef Residue res
        cdef res_iterator it
        it = self.thisptr.ResStart()

        while it != self.thisptr.ResEnd():
            res = Residue()
            res.thisptr[0] = deref(it)
            yield res
            incr(it)
        
    def mol_iter(self):
        cdef Molecule mol
        cdef mol_iterator it
        it = self.thisptr.MolStart()

        while it != self.thisptr.MolEnd():
            mol = Molecule()
            mol.thisptr[0] = deref(it)
            yield mol
            incr(it)

    def _set_parm_name(self, string title, FileName filename):
        self.thisptr.SetParmName(title, filename.thisptr[0])

    def set_reference_frame(self, Frame frame):
        """set reference frame for distance-based atommask selection

        Examples
            top.set_reference_frame(frame)
            top.select(":3 < :5.0") # select all atoms within 5.0 A to residue 3
        """
        self.thisptr.SetReferenceCoords(frame.thisptr[0])

    def file_path(self):
        return self.thisptr.c_str()

    def trunc_res_atom_name(self, id_or_mask):
        """return str or list iterator of str with format like "TYR_3@CA"
        Parameters
        ---------
        id_or_mask : int or str (no default)
            int : return single str for atom with specific ID
            str : return a list iterator of str for atoms with given mask

        """
        cdef int index 
        cdef AtomMask atm
        cdef list namelist = []

        if is_int(id_or_mask):
            index = <int> id_or_mask
            return self.thisptr.TruncResAtomName(index)
        elif isinstance(id_or_mask, string_types):
            atm = self(id_or_mask)
            for index in atm.selected_indices():
                namelist.append(self.trunc_res_atom_name(index))
            return namelist

    def atom_mask_name(self, int atom):
        return self.thisptr.AtomMaskName(atom)

    def trunc_resname_num(self, int res):
        return self.thisptr.TruncResNameNum(res)

    def find_atom_in_residue(self, int res, atname):
        cdef NameType _atomnametype
        if isinstance(atname, string_types):
            _atomnametype = NameType(atname)
        elif isinstance(atname, NameType):
            _atomnametype = <NameType> atname 
        return self.thisptr.FindAtomInResidue(res, _atomnametype.thisptr[0])
    
    @property
    def atomlist(self):
        return list(self.atoms)

    @property
    def residuelist(self):
        return list(self.residues)

    @property
    def moleculelist(self):
        return list(self.mols)

    def summary(self):
        set_world_silent(False)
        self.thisptr.Summary()
        set_world_silent(True)

    def atom_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintAtomInfo(maskString)
        set_world_silent(True)

    def bond_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintBondInfo(maskString)
        set_world_silent(True)
    
    def angle_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintAngleInfo(maskString)
        set_world_silent(True)

    def dihedral_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintDihedralInfo(maskString)
        set_world_silent(True)

    def molecule_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintMoleculeInfo(maskString)
        set_world_silent(True)

    def residue_info(self, maskString="*"):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintResidueInfo(maskString)
        set_world_silent(True)

    def charge_mass_info(self, maskString, int idtype):
        set_world_silent(False)
        maskString = maskString.encode()
        self.thisptr.PrintChargeMassInfo(maskString, idtype)
        set_world_silent(True)

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

    def set_offset(self, double x):
        self.thisptr.SetOffset(x)

    def set_ipol(self, int id):
        self.thisptr.SetIpol(id)

    @property
    def filename(self):
        # I want to keep _original_filename so don't need to
        # change other codes
        return self._original_filename

    @property
    def _original_filename(self):
        cdef FileName filename = FileName()
        filename.thisptr[0] = self.thisptr.OriginalFilename()
        return filename.__str__()

    property parm_index:
        def __get__(self):
            return self.thisptr.Pindex()

    property _p_index:
        # shortcut of parm_index
        def __get__(self):
            return self.thisptr.Pindex()

    property n_atoms:
        def __get__(self):
            return self.thisptr.Natom()

    property n_res:
        # shortcut
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

    property _parm_name:
        def __get__(self):
            return self.thisptr.ParmName()
        def __set__(self, name):
            # TODO : check
            self.thisptr.SetParmName(name, FileName().thisptr[0])

    property gb_radii:
        def __get__(self):
            return self.thisptr.GBradiiSet()

    def set_integer_mask(self, AtomMask atm, Frame frame=Frame()):
        if frame.is_empty():
            return self.thisptr.SetupIntegerMask(atm.thisptr[0])
        else:
            return self.thisptr.SetupIntegerMask(atm.thisptr[0], frame.thisptr[0])

    def _scale_dihedral_k(self, double value):
        self.thisptr.ScaleDihedralK(value)

    property box:
        def __get__(self):
            cdef Box box = Box()
            box.thisptr[0] = self.thisptr.ParmBox()
            return box
        def __set__(self, box_or_array):
            cdef Box _box
            if  isinstance(box_or_array, Box):
                _box = box_or_array
            else:
                # try to create box
                _box = Box(box_or_array)
            self.thisptr.SetParmBox(_box.thisptr[0])

    def has_box(self):
        return self.box.has_box()

    def _partial_modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.partialModifyStateByMask(m.thisptr[0]))
        return top

    def _modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.modifyStateByMask(m.thisptr[0]))
        return top

    def _modify_by_map(self, vector[int] m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.ModifyByMap(m))
        return top

    def strip_atoms(Topology self, mask, copy=False):
        """strip atoms with given mask"""
        cdef AtomMask atm
        cdef Topology new_top

        atm = AtomMask(mask)
        self.set_integer_mask(atm)
        if atm.n_atoms == 0:
            raise ValueError("number of stripped atoms must be > 1")
        atm.invert_mask()
        new_top = self._modify_state_by_mask(atm)

        if copy:
            return new_top
        else:
            self.thisptr[0] = new_top.thisptr[0]

    def is_empty(self):
        s = self.file_path()
        return s == ""

    def atom_indices(self, mask, *args, **kwd):
        """return atom indices with given mask
        To be the same as cpptraj/Ambertools: we mask indexing starts from 1
        but the return list/array use 0

        Parameters
        ---------
        mask : str
            Atom mask

        Returns
        ------
        indices : Python array
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        return atm.indices

    @property
    def atom_names(self):
        """return unique atom name in Topology
        """
        s = set()
        for atom in self.atom_iter():
            s.add(atom.name)
        return s

    @property
    def residue_names(self):
        """return unique residue names in Topology
        """
        s = set()
        for residue in self.residue_iter():
            s.add(residue.name)
        return s

    def get_parm_coord_info(self):
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

    @property
    def mass(self):
        """return python array of atom masses"""
        cdef pyarray marray = pyarray('d', [])
        cdef Atom atom

        for atom in self.atom_iter():
            marray.append(atom.mass)
        return marray

    def indices_bonded_to(self, atom_name):
        """return indices of the number of atoms that each atom bonds to
        Parameters
        ----------
        atom_name : name of the atom
        """
        cdef pyarray arr0 = pyarray('i', [])
        cdef int i, count=0

        # convert to lower case
        atom_name = atom_name.upper()

        for atom in self:
            bond_indices = atom.bonded_indices()
            count = 0
            for i in bond_indices:
                if self[i].name.startswith(atom_name):
                    count += 1
            arr0.append(count)
        return arr0

    def add_bonds(self, cython.integral [:, ::1] indices):
        """add bond for pairs of atoms. 

        Parameters
        ---------
        bond_indices : 2D array_like (must have buffer interface)
            shape=(n_atoms, 2)
        """
        cdef int i
        cdef int j, k

        for i in range(indices.shape[0]):
            j, k = indices[i, :]
            self.thisptr.AddBond(j, k)

    def add_angles(self, cython.integral [:, ::1] indices):
        """add angle for a group of 3 atoms. 

        Parameters
        ---------
        indices : 2D array_like (must have buffer interface),
            shape=(n_atoms, 3)
        """
        cdef int i
        cdef int j, k, n

        for i in range(indices.shape[0]):
            j, k, n = indices[i, :]
            self.thisptr.AddAngle(j, k, n)

    def add_dihedrals(self, cython.integral [:, ::1] indices):
        """add dihedral for a group of 4 atoms. 

        Parameters
        ---------
        indices : 2D array_like (must have buffer interface),
            shape=(n_atoms, 3)
        """
        cdef int i
        cdef int j, k, n, m

        for i in range(indices.shape[0]):
            j, k, n, m = indices[i, :]
            self.thisptr.AddDihedral(j, k, n, m)

    @property
    def bonds(self):
        """return bond iterator"""
        # both noh and with-h bonds
        cdef BondArray bondarray, bondarray_h
        cdef BondType btype = BondType()

        bondarray = self.thisptr.Bonds()
        bondarray_h = self.thisptr.BondsH()
        bondarray.insert(bondarray.end(), bondarray_h.begin(), bondarray_h.end())

        for btype.thisptr[0] in bondarray:
            yield btype

    @property
    def angles(self):
        """return bond iterator"""
        cdef AngleArray anglearray, anglearray_h
        cdef AngleType atype = AngleType()

        anglearray = self.thisptr.Angles()
        anglearray_h = self.thisptr.AnglesH()
        anglearray.insert(anglearray.end(), anglearray_h.begin(), anglearray_h.end())

        for atype.thisptr[0] in anglearray:
            yield atype

    @property
    def dihedrals(self):
        """return dihedral iterator"""
        cdef DihedralArray dharr, dharr_h
        cdef DihedralType dhtype = DihedralType()

        dharr = self.thisptr.Dihedrals()
        dharr_h = self.thisptr.DihedralsH()
        dharr.insert(dharr.end(), dharr_h.begin(), dharr_h.end())

        for dhtype.thisptr[0] in dharr:
            yield dhtype

    @property
    def bond_indices(self):
        """return an iterator of bond indices"""
        for b in self.bonds:
            yield b.indices

    @property
    def angle_indices(self):
        """return an iterator of bond indices"""
        for b in self.angles:
            yield b.indices

    @property
    def dihedral_indices(self):
        """return an iterator of bond indices"""
        for b in self.dihedrals:
            yield b.indices

    @property
    def _bonds_ndarray(self):
        _, np = _import_numpy()
        return np.asarray([b for b in self.bond_indices], dtype=np.int64)

    @property
    def _angles_ndarray(self):
        _, np = _import_numpy()
        return np.asarray([b for b in self.angle_indices], dtype=np.int64)

    @property
    def _dihedrals_ndarray(self):
        _, np = _import_numpy()
        return np.asarray([b for b in self.dihedral_indices], dtype=np.int64)

    def vdw_radii(self):
        cdef int n_atoms = self.n_atoms
        cdef int i
        cdef pyarray arr = pyarray_master.clone(pyarray('d', []), 
                           n_atoms, zero=True)
        cdef double[:] d_view = arr

        for i in range(n_atoms):
            d_view[i] = self.thisptr.GetVDWradius(i)
        return arr
