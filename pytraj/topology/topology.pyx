# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function, absolute_import
import os
import sys
cimport cython
from cython.operator cimport dereference as deref, preincrement as incr
from libcpp.string cimport string
from cpython.array cimport array as pyarray
from cpython cimport array as pyarray_master
from pytraj.c_options import set_world_silent  # turn on and off cpptraj's stdout

from collections import namedtuple
import numpy as np

from pytraj.c_dict import get_key, AtomicElementDict
from pytraj.utils.check_and_assert import is_int, is_array
from pytraj.compat import set
from pytraj.externals.six import PY2, PY3, string_types
from pytraj.externals.six.moves import range
from pytraj.utils.check_and_assert import is_int
from pytraj.c_dict import ParmFormatDict
from pytraj.utils.convert import array_to_cpptraj_atommask

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str
else:
    string_types = basestring

__all__ = ['Topology', 'ParmFile', 'SimplifiedTopology', 'SimplifiedAtom', 'SimplifiedResidue']


class SimplifiedAtom(
    namedtuple(
        'SimplifiedAtom',
        'name type element charge mass index atomic_number resname resid molnum')):
    __slots__ = ()

    def __str__(self):
        return 'SimplifiedAtom(name={}, type={}, element={}, atomic_number={}, index={}, resname={}, resid={}, molnum={})'.format(
            self.name, self.type, self.element, self.atomic_number, self.index, self.resname, self.resid, self.molnum)

    def __repr__(self):
        return str(self)


class SimplifiedResidue(
    namedtuple(
        'SimplifiedResidue',
        'name index first last')):
    __slots__ = ()

    def __str__(self):
        return 'SimplifiedResidue(name={}, index={}, atom_range={}-{})'.format(
            self.name, self.index, self.first, self.last)

    def __repr__(self):
        return str(self)


class SimplifiedTopology(namedtuple('SimplifiedTopology', 'atoms residues')):
    '''a lightweight Topology for fast iterating and convenient accessing atom, residue

    Notes
    -----
    cpptraj does not understand this class (use :class:Topology)
    '''
    __slots__ = ()

    def __str__(self):
        return 'SimplifiedTopology({} atoms, {} residues)'.format(
            len(self.atoms), len(self.residues))

    def __repr__(self):
        return str(self)


cdef class Topology:
    def __cinit__(self, *args):
        """
        args = filename or Topology instance
        """
        cdef Topology tp
        cdef string filename
        self.thisptr = new _Topology()

        # I dont make default _own_memory (True) in __cinit__ since
        # when passing something like top = Topology(filename), Python/Cython
        # would think "filename" is value of _own_memory
        self._own_memory = True

        if not args:
            # print "there is no args" # for debug
            # make empty Topology instance
            pass
        else:
            if len(args) == 1:
                if isinstance(args[0], Topology):
                    tp = args[0]
                    self.thisptr[0] = tp.thisptr[0]
                else:
                    raise ValueError()
            else:
                raise ValueError()

    def __dealloc__(self):
        if self._own_memory and self.thisptr:
            del self.thisptr

    def __str__(self):
        box = self.box
        if box.has_box():
            box_txt = "PBC with box type = %s" % box.type
        else:
            box_txt = "non-PBC"

        tmp = "<%s: %s atoms, %s residues, %s mols, %s>" % (
            self.__class__.__name__,
            self.n_atoms,
            self.n_residues,
            self.n_mols,
            box_txt)
        return tmp

    def add_atom(self, Atom atom, Residue residue):
        self.thisptr.AddTopAtom(atom.thisptr[0], residue.thisptr[0])

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.n_atoms

    def __mul__(self, int n_times):
        cdef int i
        t = self.copy()

        for i in range(n_times - 1):
            t.join(self.copy())
        return t

    def __add__(self, Topology other):
        '''order matters

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj.testing import get_fn
        >>> tn0 = get_fn('tz2')[1]
        >>> tn1 = get_fn('ala3')[1]
        >>> top0 = pt.load_topology(tn0)
        >>> top0
        <Topology: 5293 atoms, 1704 residues, 1692 mols, PBC with box type = ortho>
        >>> top1 = pt.load_topology(tn1)
        >>> top1
        <Topology: 34 atoms, 3 residues, 1 mols, non-PBC>
        >>> top0 + top1
        <Topology: 5327 atoms, 1707 residues, 1693 mols, PBC with box type = ortho>
        '''
        new_top = self.copy()
        new_top.join(other)
        return new_top

    def __iadd__(self, Topology other):
        self.join(other)
        return self

    def load(self, string filename):
        """loading Topology from filename. This is for internal use. Should ``pytraj.load_topology``

        """
        del self.thisptr
        self = Topology(filename)

    def set_distance_mask_reference(self, Frame frame):
        self.thisptr.SetDistMaskRef(frame.thisptr[0])

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
        """return an Atom, a list of Atom or a new Topology

        Returns
        -------
        an Atom if idx is an integer
        a new Topology if idx a string mask
        other cases: a list of Atoms

        Examples
        --------
        In [31]: top[0]
        Out[31]: <N-atom, resid=0, n_bonds=4>
        """

        cdef Atom atom
        cdef AtomMask atm
        cdef Residue res
        cdef int i
        cdef list alist = []

        if is_int(idx):
            # need to explicitly cast to int
            i = <int > idx
            atom = Atom()
            if i >= 0:
                atom.thisptr[0] = self.thisptr.index_opr(i)
            else:
                # negative indexing
                atom.thisptr[0] = self.thisptr.index_opr(self.n_atoms + i)
            return atom
        elif isinstance(idx, string_types):
            # return atom object iterator with given mask
            # self(idx) return AtomMask object
            mask = idx
        elif isinstance(idx, AtomMask):
            atm = <AtomMask > idx
            # return atom object iterator with given mask
            # self(idx) return AtomMask object
            mask = array_to_cpptraj_atommask(idx.indices)
        elif isinstance(idx, slice):
            # does not have memory efficiency with large Topology
            # (since we convert to atom list first)
            start, stop, step = idx.indices(self.n_atoms)
            mask = array_to_cpptraj_atommask(range(start, stop, step))
        elif isinstance(idx, (list, tuple, range)) or is_array(idx):
            mask = array_to_cpptraj_atommask(idx)
        elif isinstance(idx, Residue):
            mask = array_to_cpptraj_atommask(
                range(idx.first_atom_index, idx.last_atom_index))
            return self._get_new_from_mask(mask)
        elif isinstance(idx, Molecule):
            mol = idx
            mask = array_to_cpptraj_atommask(range(mol.begin_atom, mol.end_atom))
        else:
            raise NotImplementedError("")

        return self._get_new_from_mask(mask)

    def __call__(self, mask, *args, **kwd):
        """intended to use with Frame indexing: atm = top('@CA') (for internal use)
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        return atm

    def __iter__(self):
        return self.atoms

    def _iter_mut(self):
        cdef Atom atom
        cdef int n_atoms = self.n_atoms
        cdef int idx = 0

        for idx in range(n_atoms):
            atom = Atom()
            atom.thisptr = &(self.thisptr.GetAtomView(idx))
            atom.own_memory = False
            atom.index = idx
            # do not call python object here to avoid overhead
            #atom.residue = self._residue_light(atom.resid)
            yield atom

    def select(self, mask):
        """return atom indices

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2()
        >>> atm = traj.top.select("@CA")
        >>> atm
        array([  4,  15,  39, ..., 159, 173, 197])
        >>> pt.rmsd(traj, mask=atm)
        array([  1.94667955e-07,   2.54596866e+00,   4.22333034e+00, ...,
                 4.97189564e+00,   5.53947712e+00,   4.83201237e+00])

        Notes
        -----
        support openmp for distance-based atommask selction
        """
        return self(mask).indices

    property atoms:
        def __get__(self):
            cdef Atom atom
            cdef atom_iterator it
            cdef int idx = 0

            it = self.thisptr.begin()
            while it != self.thisptr.end():
                atom = Atom()
                atom.thisptr[0] = deref(it)
                atom.index = idx
                atom.resname = self.thisptr.Res(atom.resid).c_str().strip()
                yield atom
                idx += 1
                incr(it)

    def simplify(self):
        '''return a (immutable) light version of Topology for fast iterating. (experiment)

        No writing capabibility (you should use ParmEd for Topology editing)

        Examples
        --------
        >>> import pytraj as pt
        >>> top = pt.load_topology('data/tz2.parm7')

        >>> simp_top = top.simplify()
        >>> atom = simp_top.atoms[0]
        >>> atom.resname
        'SER'

        >>> res = simp_top.residues[0]
        >>> # get all atoms for 1st residue
        >>> atoms = simp_top.atoms[res.first:res.last]
        '''
        cdef _Atom atom
        cdef atom_iterator ait
        cdef res_iterator rit
        cdef int idx = 0
        cdef _Residue res
        cdef list atoms, residues

        atoms = []
        residues = []

        # get atoms
        ait = self.thisptr.begin()
        while ait != self.thisptr.end():
            atom = deref(ait)
            res = self.thisptr.Res(atom.ResNum())
            atoms.append(SimplifiedAtom(name=atom.c_str().strip(),
                                        type=atom.Type().Truncated(),
                                        element=get_key(
                                            atom.Element(), AtomicElementDict).lower(),
                                        charge=atom.Charge(),
                                        mass=atom.Mass(),
                                        index=idx,
                                        atomic_number=atom.AtomicNumber(),
                                        resname=res.c_str().strip(),
                                        resid=res.OriginalResNum() - 1,
                                        molnum=atom.MolNum()))
            idx += 1
            incr(ait)

        # get residues
        rit = self.thisptr.ResStart()
        idx = 0
        while rit != self.thisptr.ResEnd():
            res = deref(rit)
            residues.append(SimplifiedResidue(name=res.c_str().strip(),
                                              index=idx,
                                              first=res.FirstAtom(),
                                              last=res.LastAtom()))
            idx += 1
            incr(rit)
        return SimplifiedTopology(atoms=atoms, residues=residues)

    property residues:
        def __get__(self):
            cdef Residue res
            cdef res_iterator it
            it = self.thisptr.ResStart()

            while it != self.thisptr.ResEnd():
                res = Residue()
                res.thisptr[0] = deref(it)
                yield res
                incr(it)

    property mols:
        def __get__(self):
            cdef Molecule mol
            cdef mol_iterator it
            it = self.thisptr.MolStart()

            while it != self.thisptr.MolEnd():
                mol = Molecule()
                mol.thisptr[0] = deref(it)
                yield mol
                incr(it)

    property atomlist:
        '''return a copy of atoms. If the Topology is large, this method calling
        is every expensive. Make sure to call once and save it.
        '''

        def __get__(self):
            return list(self.atoms)

    property residuelist:
        """return a copy of residues
        """
        def __get__(self):
            return list(self.residues)

    property moleculelist:
        """return a copy of molecules. (not much information)
        """
        def __get__(self):
            return list(self.mols)

    def summary(self):
        """basic info. This information only appears in Ipython or Python shell. 
        It does not appear in Jupyter notebook (due to C++ stdout)
        """
        set_world_silent(False)
        self.thisptr.Summary()
        set_world_silent(True)

    def start_new_mol(self):
        self.thisptr.StartNewMol()

    property filename:
        """return original filename. This is for testing purpose.
        """
        def __get__(self):
            # I want to keep _original_filename so don't need to
            # change other codes
            import os
            return os.path.abspath(self._original_filename)

    property _original_filename:
        def __get__(self):
            '''NOTE: do not delete me
            '''
            cdef FileName filename = FileName()
            filename.thisptr[0] = self.thisptr.OriginalFilename()
            return filename.__str__()

    property n_atoms:
        def __get__(self):
            return self.thisptr.Natom()

    property n_residues:
        def __get__(self):
            return self.thisptr.Nres()

    property n_mols:
        def __get__(self):
            return self.thisptr.Nmol()

    property n_solvents:
        def __get__(self):
            return self.thisptr.Nsolvent()

    def set_integer_mask(self, AtomMask atm, Frame frame=Frame()):
        if frame.n_atoms == 0:
            return self.thisptr.SetupIntegerMask(atm.thisptr[0])
        else:
            return self.thisptr.SetupIntegerMask(atm.thisptr[0], frame.thisptr[0])

    property box:
        def __get__(self):
            cdef Box box = Box()
            box.thisptr[0] = self.thisptr.ParmBox()
            return box

        def __set__(self, box_or_array):
            cdef Box _box
            if isinstance(box_or_array, Box):
                _box = box_or_array
            else:
                # try to create box
                _box = Box(box_or_array)
            self.thisptr.SetParmBox(_box.thisptr[0])

    def has_box(self):
        return self.box.has_box()

    def set_nobox(self):
        self.box = Box()
        return self

    def _partial_modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.partialModifyStateByMask(m.thisptr[0]))
        return top

    def _modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.thisptr.modifyStateByMask(m.thisptr[0]))
        return top

    def _get_new_from_mask(self, mask=None):
        '''

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2()
        >>> top = traj.top
        >>> top._get_new_from_mask('@CA')
        <Topology: 12 atoms, 12 residues, 12 mols, non-PBC>
        '''
        if mask is None or mask == "":
            return self
        else:
            atm = self(mask)
            return self._modify_state_by_mask(atm)

    def strip(Topology self, mask, copy=False):
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
        return self.n_atoms == 0

    def atom_indices(self, mask, *args, **kwd):
        """return atom indices with given mask
        To be the same as cpptraj/Ambertools: we mask indexing starts from 1
        but the return list/array use 0

        Parameters
        ----------
        mask : str
            Atom mask

        Returns
        ------
        indices : Python array
        """
        cdef AtomMask atm = AtomMask(mask)
        self.set_integer_mask(atm)
        return atm.indices

    property atom_names:
        """return unique atom name in Topology
        """
        def __get__(self):
            s = set()
            for atom in self.atoms:
                s.add(atom.name)
            return s

    property residue_names:
        """return unique residue names in Topology
        """
        def __get__(self):
            s = set()
            for residue in self.residues:
                s.add(residue.name)
            return s

    def join(self, Topology top):
        if top is self:
            raise ValueError('must not be your self')
        self.thisptr.AppendTop(top.thisptr[0])
        return self

    property mass:
        '''return a copy of atom masses (numpy 1D array)'''

        def __get__(self):
            """return python array of atom masses"""
            cdef Atom atom
            return np.asarray([atom.mass for atom in self.atoms])

    property charge:
        '''return a copy of atom charges (numpy 1D array)'''

        def __get__(self):
            return np.asarray([atom.charge for atom in self.atoms])

    def indices_bonded_to(self, atom_name):
        """return indices of the number of atoms that each atom bonds to

        Parameters
        ----------
        atom_name : name of the atom
        """
        cdef pyarray arr0 = pyarray('i', [])
        cdef int i, count = 0

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

    def add_bonds(self, cython.integral[:, ::1] indices):
        """add bond for pairs of atoms.

        Parameters
        ----------
        bond_indices : 2D array_like (must have buffer interface)
            shape=(n_atoms, 2)
        """
        cdef int i
        cdef int j, k

        for i in range(indices.shape[0]):
            j, k = indices[i, :]
            self.thisptr.AddBond(j, k)

    def add_angles(self, cython.integral[:, ::1] indices):
        """add angle for a group of 3 atoms.

        Parameters
        ----------
        indices : 2D array_like (must have buffer interface),
            shape=(n_atoms, 3)
        """
        cdef int i
        cdef int j, k, n

        for i in range(indices.shape[0]):
            j, k, n = indices[i, :]
            self.thisptr.AddAngle(j, k, n)

    def add_dihedrals(self, cython.integral[:, ::1] indices):
        """add dihedral for a group of 4 atoms.

        Parameters
        ----------
        indices : 2D array_like (must have buffer interface),
            shape=(n_atoms, 3)
        """
        cdef int i
        cdef int j, k, n, m

        for i in range(indices.shape[0]):
            j, k, n, m = indices[i, :]
            self.thisptr.AddDihedral(j, k, n, m)

    property bonds:
        def __get__(self):
            """return bond iterator"""
            # both noh and with-h bonds
            cdef BondArray bondarray, bondarray_h
            cdef BondType btype = BondType()

            bondarray = self.thisptr.Bonds()
            bondarray_h = self.thisptr.BondsH()
            bondarray.insert(bondarray.end(), bondarray_h.begin(), bondarray_h.end())

            for btype.thisptr[0] in bondarray:
                yield btype

    property angles:
        def __get__(self):
            """return bond iterator"""
            cdef AngleArray anglearray, anglearray_h
            cdef AngleType atype = AngleType()

            anglearray = self.thisptr.Angles()
            anglearray_h = self.thisptr.AnglesH()
            anglearray.insert(anglearray.end(), anglearray_h.begin(), anglearray_h.end())

            for atype.thisptr[0] in anglearray:
                yield atype

    property dihedrals:
        def __get__(self):
            """return dihedral iterator"""
            cdef DihedralArray dharr, dharr_h
            cdef DihedralType dhtype = DihedralType()

            dharr = self.thisptr.Dihedrals()
            dharr_h = self.thisptr.DihedralsH()
            dharr.insert(dharr.end(), dharr_h.begin(), dharr_h.end())

            for dhtype.thisptr[0] in dharr:
                yield dhtype

    property bond_indices:
        def __get__(self):
            return np.asarray([b.indices for b in self.bonds], dtype=np.int64)

    property angle_indices:
        def __get__(self):
            return np.asarray([b.indices for b in self.angles], dtype=np.int64)

    property dihedral_indices:
        def __get__(self):
            return np.asarray([b.indices for b in self.dihedrals], dtype=np.int64)

    def __getstate__(self):
        return self.to_dict()

    def __setstate__(self, dict_data):
        d = dict_data

        # always start molnum at 0.
        MOLNUM = 0

        for idx, (aname, atype, charge, mass, resid, resname, mol_number) in enumerate(zip(d['atom_name'], d[
                'atom_type'], d['atom_charge'], d['atom_mass'], d['resid'], d['resname'], d['mol_number'])):
            atom = Atom(name=aname, type=atype, charge=charge, mass=mass, resid=resid)
            atom.set_mol(mol_number)
            residue = Residue(resname, resid)
            if idx == 0:
                self.start_new_mol()
            if mol_number > MOLNUM:
                self.start_new_mol()
                MOLNUM += 1
            self.add_atom(atom, residue)

        # add box
        box = Box(d['box'])
        self.box = box

        self.add_bonds(d['bond_index'])
        self.add_dihedrals(d['dihedral_index'])

    @classmethod
    def from_dict(cls, dict_data):
        """internal use for serialize Topology
        """
        new_top = Topology()
        new_top.__setstate__(dict_data)
        return new_top

    def to_dict(self):
        '''convert Topology to Python dict
        '''
        cdef:
            int n_atoms = self.n_atoms
            int idx
            Atom atom

        d = {}

        short_resnamelist = np.asarray([res.name for res in self.residues])
        resids = []

        atomnames = []
        atomtypes = []
        atomcharges = []
        molnums = []
        resnames = []

        for idx, atom in enumerate(self.atoms):
            resids.append(atom.resid)
            atomnames.append(atom.name)
            atomtypes.append(atom.type)
            atomcharges.append(atom.charge)
            molnums.append(atom.molnum)
            resnames.append(short_resnamelist[atom.resid])

        d['atom_name'] = atomnames
        d['atom_type'] = atomtypes
        d['atom_charge'] = atomcharges
        d['atom_mass'] = self.mass
        d['resname'] = resnames
        d['resid'] = resids
        d['bond_index'] = self.bond_indices
        d['dihedral_index'] = self.dihedral_indices
        d['mol_number'] = molnums
        d['box'] = self.box.values

        return d

    def to_dataframe(self):
        """convert to pandas' DataFrame. (experiment)
        """
        import pandas as pd
        cdef:
            int n_atoms = self.n_atoms
            int idx
            Atom atom

        if pd:
            labels = ['resid', 'resname', 'atomname', 'atomic_number', 'mass']
            mass_arr = np.array(self.mass)
            resid_arr = np.empty(n_atoms, dtype='i')
            resname_arr = np.empty(n_atoms, dtype='U4')
            atomname_arr = np.empty(n_atoms, 'U4')
            atomicnumber_arr = np.empty(n_atoms, dtype='i4')

            for idx, atom in enumerate(self.atoms):
                # TODO: make faster?
                resid_arr[idx] = atom.resid
                resname_arr[idx] = self.residuelist[atom.resid].name
                atomname_arr[idx] = atom.name
                atomicnumber_arr[idx] = atom.atomic_number

            arr = np.vstack((resid_arr, resname_arr, atomname_arr,
                             atomicnumber_arr, mass_arr)).T
            return pd.DataFrame(arr, columns=labels)
        else:
            raise ValueError("must have pandas")

    def to_parmed(self):
        """try to load to ParmEd's Structure
        """
        import parmed as pmd
        from pytraj.utils import tempfolder

        with tempfolder():
            self.save("tmp.prmtop", overwrite=True)
            return pmd.load_file("tmp.prmtop")

    property _total_charge:
        def __get__(self):
            return sum([atom.charge for atom in self.atoms])

    def save(self, filename=None, format='AMBERPARM'):
        """save to given file format (parm7, psf, ...)
        """
        parm = ParmFile()
        parm.writeparm(filename=filename, top=self, format=format)

    def set_solvent(self, mask):
        '''set ``mask`` as solvent
        '''
        mask = mask.encode()
        self.thisptr.SetSolvent(mask)

    def residue(self, int idx, bint atom=False):
        '''
        if atom is True, get full list of atoms for idx-th residue. This will be very slow
        if atom is False, get ()
        '''
        cdef Residue res = Residue()
        res.thisptr[0] = self.thisptr.Res(idx)
        return res

    def atom(self, int idx):
        '''return an Atom based on idx. Update this Atom will update Topology.
        Make this method private for now.
        '''
        cdef Atom atom = Atom()
        atom.own_memory = False
        atom.thisptr = &self.thisptr.GetAtomView(idx)
        atom.resname = self.thisptr.Res(atom.resid).c_str().strip()
        return atom


cdef class ParmFile:
    def __cinit__(self):
        self.thisptr = new _ParmFile()

    def __dealloc__(self):
        del self.thisptr

    def readparm(self, filename="", top=Topology(), option=''):
        """readparm(Topology top=Topology(), string filename="", "*args)
        Return : None (update `top`)

        top : Topology instance
        filename : str, output filename
        arglist : ArgList instance, optional

        """
        cdef ArgList arglist
        cdef debug = 0
        cdef Topology _top = <Topology > top

        filename = filename.encode()

        if not option:
            self.thisptr.ReadTopology(_top.thisptr[0], filename, debug)
        else:
            if isinstance(option, string_types):
                arglist = ArgList(option)
            else:
                arglist = <ArgList > option
            self.thisptr.ReadTopology(
                _top.thisptr[0], filename, arglist.thisptr[0], debug)

    def writeparm(self, Topology top=Topology(), filename="default.top",
                  ArgList arglist=ArgList(), format=""):
        cdef int debug = 0
        cdef int err
        # change `for` to upper
        cdef ParmFormatType parmtype
        filename = filename.encode()

        if format == "":
            parmtype = UNKNOWN_PARM
        else:
            try:
                format = format.upper()
                parmtype = ParmFormatDict[format]
            except:
                raise ValueError("supported keywords: ", self.formats)

        if top.is_empty():
            raise ValueError("empty topology")

        err = self.thisptr.WriteTopology(
            top.thisptr[0], filename, arglist.thisptr[0], parmtype, debug)
        if err == 1:
            raise ValueError("Not supported or failed to write")

    def filename(self):
        cdef FileName filename = FileName()
        filename.thisptr[0] = self.thisptr.ParmFilename()
        return os.path.abspath(filename)
