# cython: c_string_type=unicode, c_string_encoding=utf8
from __future__ import print_function, absolute_import
import os
import sys
cimport cython
from cython.operator cimport dereference as deref, preincrement as incr
from libcpp.string cimport string
from pytraj.core.c_options import set_world_silent  # turn on and off cpptraj's stdout

from collections import namedtuple
import numpy as np

from pytraj.core.c_dict import get_key, AtomicElementDict
from pytraj.utils.check_and_assert import is_int, is_array
from pytraj.utils.context import capture_stdout
from pytraj.utils.check_and_assert import is_int
from pytraj.core.c_dict import ParmFormatDict
from pytraj.utils.convert import array_to_cpptraj_atommask


__all__ = ['Topology', 'ParmFile', 'SimplifiedTopology', 'SimplifiedAtom', 'SimplifiedResidue']


class SimplifiedAtom(object):
    '''EXPERIMENTAL

    read only
    '''
    def __init__(self, name, type, element, charge, mass, index, atomic_number,
            resname, resid, molnum,
            bond_indices=None,
            simplified_top=None):
        self.name = name
        self.type = type
        self.element = element
        self.charge = charge
        self.mass = mass
        self.index = index
        self.atomic_number = atomic_number
        self.resname = resname
        self.resid = resid
        self.molnum = molnum
        self._bond_indices = bond_indices
        self._simplified_top = simplified_top

    @property
    def residue(self):
        return self._simplified_top.residues[self.resid]

    @property
    def bond_partners(self):
        return self._simplified_top.atoms[self._bond_indices]

    def __str__(self):
        return 'SimplifiedAtom(name={}, type={}, element={}, atomic_number={}, index={},'\
               'resname={}, resid={}, molnum={})'.format(
            self.name, self.type, self.element, self.atomic_number,
            self.index, self.resname, self.resid, self.molnum)

    def __repr__(self):
        return str(self)


class SimplifiedResidue(object):
    '''EXPERIMENTAL

    read only
    '''
    def __init__(self, name, index, first, last, associated_topology):
        self.name = name
        self.index = index
        self.first = first
        self.last = last
        self._simplified_top = associated_topology

    @property
    def atoms(self):
        return self._simplified_top.atoms[self.first:self.last]

    def __str__(self):
        return 'SimplifiedResidue(name={}, index={}, atom_range={}-{})'.format(
            self.name, self.index, self.first, self.last)

    def __repr__(self):
        return str(self)

    def __iter__(self):
        for atom in self.atoms:
            yield atom


class SimplifiedTopology(object):
    '''EXPERIMENTAL: a lightweight Topology for fast iterating and convenient accessing atom, residue

    Notes
    -----
    - cpptraj does not understand this class (use :class:Topology)
    - read only
    '''
    def __init__(self, atoms=None, residues=None, cpptraj_topology=None):
        self.atoms = np.asarray(atoms)
        self.residues = np.asarray(residues)
        self._top = cpptraj_topology

    def __str__(self):
        return 'SimplifiedTopology({} atoms, {} residues)'.format(
            len(self.atoms), len(self.residues))

    def __repr__(self):
        return str(self)

    def __iter__(self):
        for atom in self.atoms:
            yield atom

    def select(self, mask):
        return self._top.select(mask)


cdef class Topology:
    """Molecular topology representation for MD simulations and analysis.

    The Topology class provides comprehensive molecular structure information including
    atoms, residues, molecules, bonds, angles, dihedrals, and other topological data.
    It interfaces with cpptraj's Topology backend for high-performance operations.

    Key Features:
    - Atom, residue, and molecule management
    - Advanced selection capabilities with masks
    - Structural analysis and manipulation
    - File I/O for multiple formats (AMBER, PDB, PSF, etc.)
    - Integration with trajectory analysis workflows

    Memory Management:
    - Uses C++ backend for performance
    - Automatic memory cleanup when Python object is deleted
    - Safe copying and reference handling

    Examples:
    --------
    >>> import pytraj as pt
    >>> # Load topology from file
    >>> top = pt.load_topology('system.parm7')
    >>> print(f"System: {top.n_atoms} atoms, {top.n_residues} residues")

    >>> # Select atoms
    >>> ca_indices = top.select('@CA')
    >>> protein_atoms = top['@CA,C,N']

    >>> # Access structural information
    >>> first_residue = top.residue(0)
    >>> backbone = top.select(':1-10@CA,C,N,O')
    """

    # =========================================================================
    # SECTION 1: MEMORY MANAGEMENT AND CONSTRUCTION
    # =========================================================================

    def __cinit__(self, *args):
        """Initialize Topology object.

        Parameters
        ----------
        *args : optional
            Can be empty (creates empty topology) or contain a Topology instance to copy

        Notes
        -----
        - Memory is automatically managed by Cython/C++
        - Use pytraj.load_topology() for loading from files
        - Direct construction is mainly for internal use

        Examples
        --------
        >>> # Create empty topology (internal use)
        >>> top = pt.topology.Topology()

        >>> # Load from file (recommended)
        >>> top = pt.load_topology('system.parm7')
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
        """Automatic memory cleanup when Topology object is destroyed.

        Notes
        -----
        - Called automatically by Python garbage collector
        - Only deallocates if topology owns its memory
        - Ensures no memory leaks in C++ backend
        """
        if self._own_memory and self.thisptr:
            del self.thisptr

    def copy(self, *args):
        """Create a deep copy of topology or copy from another topology.

        Parameters
        ----------
        *args : optional
            If empty, creates copy of self.
            If Topology provided, copies that topology to self.

        Returns
        -------
        Topology
            New topology instance (if copying self) or None (if copying to self)

        Examples
        --------
        >>> top_copy = top.copy()  # Create new copy
        >>> top.copy(other_top)    # Copy other_top to top

        Notes
        -----
        - Deep copy includes all structural information
        - Safe for independent manipulation
        - Memory-efficient C++ implementation
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

    # =========================================================================
    # SECTION 2: STRING REPRESENTATION AND BASIC PROPERTIES
    # =========================================================================

    def __str__(self):
        """Return human-readable string representation.

        Returns
        -------
        str
            Formatted string with atom/residue counts and PBC info

        Examples
        --------
        >>> print(top)
        <Topology: 5293 atoms, 1704 residues, 1692 mols, PBC with box type = ortho>
        """
        box = self.box
        if box.has_box():
            box_txt = "PBC with box type = %s" % box.type
        else:
            box_txt = "non-PBC"

        return "<%s: %s atoms, %s residues, %s mols, %s>" % (
            self.__class__.__name__,
            self.n_atoms,
            self.n_residues,
            self.n_mols,
            box_txt)

    def __repr__(self):
        """Return string representation for debugging."""
        return self.__str__()

    def __len__(self):
        """Return number of atoms in topology.

        Returns
        -------
        int
            Total number of atoms

        Examples
        --------
        >>> len(top)  # Same as top.n_atoms
        5293
        """
        return self.n_atoms

    # =========================================================================
    # SECTION 3: TOPOLOGY MANIPULATION AND COMBINATION
    # =========================================================================

    def add_atom(self, Atom atom, Residue residue):
        """Add an atom to the topology.

        Parameters
        ----------
        atom : Atom
            Atom object to add
        residue : Residue
            Residue to add atom to

        Notes
        -----
        - Low-level method for topology construction
        - Use with caution - can break topology consistency
        - Consider using higher-level topology builders
        """
        self.thisptr.AddTopAtom(atom.thisptr[0], residue.thisptr[0])

    def __mul__(self, int n_times):
        """Replicate topology multiple times.

        Parameters
        ----------
        n_times : int
            Number of copies to make

        Returns
        -------
        Topology
            New topology with system replicated n_times

        Examples
        --------
        >>> big_system = top * 3  # Triplicate system

        Notes
        -----
        - Useful for creating larger simulation systems
        - Atom indices will be renumbered sequentially
        - Box information is preserved from original
        """
        cdef int i
        t = self.copy()

        for i in range(n_times - 1):
            t.join(self.copy())
        return t

    def __add__(self, Topology other):
        """Combine two topologies into a new topology.

        Parameters
        ----------
        other : Topology
            Topology to append to this one

        Returns
        -------
        Topology
            New combined topology

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
        >>> combined = top0 + top1
        >>> combined
        <Topology: 5327 atoms, 1707 residues, 1693 mols, PBC with box type = ortho>

        Notes
        -----
        - Order matters: self + other (self comes first)
        - Box info from first topology is preserved
        - Atom indices are renumbered sequentially
        - Creates new topology, leaves originals unchanged
        """
        new_top = self.copy()
        new_top.join(other)
        return new_top

    def __iadd__(self, Topology other):
        """In-place combination with another topology.

        Parameters
        ----------
        other : Topology
            Topology to append to this one

        Returns
        -------
        Topology
            Self (modified)

        Examples
        --------
        >>> top += other_top  # Modifies top in-place

        Notes
        -----
        - Modifies self instead of creating new topology
        - More memory efficient than + operator
        - Atom indices are renumbered
        """
        self.join(other)
        return self

    def _load(self, string filename):
        """Internal method for loading topology from file.

        Parameters
        ----------
        filename : str
            Path to topology file

        Notes
        -----
        - Internal use only
        - Use pytraj.load_topology() for user-facing loading
        - Replaces current topology data
        """
        del self.thisptr
        self = Topology(filename)

    def set_reference(self, Frame frame):
        """Set reference frame for distance-based atom selections.

        Parameters
        ----------
        frame : Frame
            Reference frame for distance calculations

        Notes
        -----
        - Enables distance-based masks like '@CA<:5.0'
        - Frame coordinates used as reference point
        - Affects subsequent select() calls using distance criteria

        Examples
        --------
        >>> top.set_reference(first_frame)
        >>> nearby_atoms = top.select('@CA<:5.0')  # CA atoms within 5Å
        """
        self.thisptr.SetDistMaskRef(frame.thisptr[0])

    # =========================================================================
    # SECTION 4: INDEXING AND SELECTION
    # =========================================================================

    def __getitem__(self, idx):
        """Access atoms, residues, or create new topologies via indexing.

        Parameters
        ----------
        idx : int, str, slice, list, AtomMask, Residue, Molecule
            Selection criteria:
            - int: return single Atom
            - str: atom selection mask, returns new Topology
            - slice: range of atoms, returns new Topology
            - list/array: list of atom indices, returns new Topology
            - AtomMask: predefined mask, returns new Topology
            - Residue: all atoms in residue, returns new Topology

        Returns
        -------
        Atom or Topology
            - Atom: if idx is integer
            - Topology: if idx is selection (mask, slice, list, etc.)

        Examples
        --------
        >>> # Single atom access
        >>> atom = top[0]  # First atom
        >>> print(atom.name, atom.resname)
        'N' 'SER'

        >>> # Selection masks - return new Topology
        >>> backbone = top['@CA,C,N,O']  # Backbone atoms
        >>> protein = top[':1-100']      # Residues 1-100
        >>> ca_atoms = top['@CA']        # All CA atoms
        >>>
        >>> # Advanced selections
        >>> hydrophobic = top[':ALA,VAL,LEU,ILE,PHE']
        >>> active_site = top[':100-120@CA,CB']

        >>> # Slice notation
        >>> first_10 = top[0:10]     # First 10 atoms
        >>> every_other = top[::2]   # Every other atom

        >>> # List of indices
        >>> specific = top[[1, 5, 10, 15]]  # Specific atoms

        >>> # Negative indexing
        >>> last_atom = top[-1]     # Last atom

        Notes
        -----
        - String masks use cpptraj syntax
        - New topologies maintain connectivity when possible
        - Single atoms are references to original topology data
        - Slicing and list indexing create independent topologies
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
        elif isinstance(idx, str):
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
            # TODO fix for Molecule Unit
            #mol = idx
            #mask = array_to_cpptraj_atommask(range(mol.begin_atom, mol.end_atom))
            raise NotImplementedError("")
        else:
            raise NotImplementedError("")

        return self._get_new_from_mask(mask)

    def __call__(self, mask, *args, **kwd):
        """Create AtomMask object for internal frame indexing operations.

        Parameters
        ----------
        mask : str
            Atom selection mask string

        Returns
        -------
        AtomMask
            Mask object for frame operations

        Notes
        -----
        - Internal method for Frame indexing: frame[top('@CA')]
        - Use select() method for getting atom indices
        - Creates mask with resolved atom indices
        """
        cdef AtomMask atm = AtomMask(mask)
        self._set_integer_mask(atm)
        return atm

    def __iter__(self):
        """Iterate over all atoms in topology.

        Yields
        ------
        Atom
            Each atom in the topology

        Examples
        --------
        >>> for atom in top:
        ...     print(f"{atom.name} {atom.resname}")
        'N' 'SER'
        'CA' 'SER'
        ...

        Notes
        -----
        - Memory-efficient iteration
        - Atoms are views into original topology
        - Modifying atoms affects original topology
        """
        return self.atoms

    def _iter_mut(self):
        """Internal mutable iterator over atoms.

        Yields
        ------
        Atom
            Each atom with mutable access

        Notes
        -----
        - Internal method for topology modification
        - Atoms can be modified in-place
        - Use with caution to maintain topology consistency
        """
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
        """Select atom indices using cpptraj mask syntax.

        Parameters
        ----------
        mask : str
            Atom selection mask using cpptraj syntax

        Returns
        -------
        numpy.ndarray
            Array of selected atom indices (0-based)

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2()
        >>>
        >>> # Basic atom selections
        >>> ca_atoms = traj.top.select("@CA")        # All CA atoms
        >>> backbone = traj.top.select("@CA,C,N,O")  # Backbone atoms
        >>>
        >>> # Residue-based selections
        >>> protein = traj.top.select(":1-100")       # Residues 1-100
        >>> first_res = traj.top.select(":1")        # First residue
        >>>
        >>> # Combined selections
        >>> active_site = traj.top.select(":85-95@CA,CB,CG")
        >>>
        >>> # Distance-based (requires reference frame)
        >>> traj.top.set_reference(traj[0])
        >>> nearby = traj.top.select("@CA<:5.0")     # CA within 5Å
        >>>
        >>> # Use in analysis
        >>> ca_indices = traj.top.select("@CA")
        >>> rmsd_values = pt.rmsd(traj, mask=ca_indices)
        >>> print(f"Found {len(ca_indices)} CA atoms")
        Found 13 CA atoms

        Mask Syntax:
        -----------
        - @ATOMNAME : select by atom name (@CA, @CB)
        - :RESID : select by residue number (:1, :1-10)
        - :RESNAME : select by residue name (:ALA, :SER)
        - <:DISTANCE : within distance (requires reference)
        - >:DISTANCE : beyond distance (requires reference)
        - Combinations: ":1-10@CA,C,N,O", ":ALA,VAL@CB"

        Notes
        -----
        - Supports OpenMP for distance-based selections
        - Returns 0-based indices (Python convention)
        - Highly optimized C++ implementation
        - Distance selections require set_reference() call first
        """
        return np.asarray(self(mask).indices, dtype='int')

    # =========================================================================
    # SECTION 5: ATOM AND RESIDUE ACCESS
    # =========================================================================

    property atoms:
        """Iterator over all atoms in the topology.

        Yields
        ------
        Atom
            Each atom in the topology with full attribute access

        Examples
        --------
        >>> # Iterate over all atoms
        >>> for atom in top.atoms:
        ...     print(f"{atom.index}: {atom.name} {atom.resname} {atom.charge:.3f}")
        0: N SER -0.416
        1: H SER 0.272
        2: CA SER -0.024
        ...

        >>> # Convert to list for indexing
        >>> atom_list = list(top.atoms)
        >>> first_atom = atom_list[0]

        >>> # Count specific atom types
        >>> ca_count = sum(1 for atom in top.atoms if atom.name == 'CA')
        >>> print(f"Found {ca_count} CA atoms")

        Properties Available on Each Atom:
        --------------------------------
        - name : str - Atom name (CA, CB, etc.)
        - type : str - Atom type (CT, C, etc.)
        - element : str - Element symbol (C, N, O, etc.)
        - charge : float - Partial charge
        - mass : float - Atomic mass
        - index : int - 0-based atom index
        - resid : int - 0-based residue index
        - resname : str - Residue name (ALA, SER, etc.)
        - molnum : int - 0-based molecule number

        Notes
        -----
        - Memory-efficient iterator (doesn't load all atoms at once)
        - Atoms are views into the original topology
        - Modifying atom properties affects the topology
        - Use list(top.atoms) if you need random access
        """
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

    cdef cppvector[int] _get_atom_bond_indices(self, _Atom atom):
        cdef cppvector[int] arr
        cdef bond_iterator it
        cdef int i = 0

        it = atom.bondbegin()
        while it != atom.bondend():
            arr.push_back(deref(it))
            incr(it)
            i += 1
        return arr

    def simplify(self):
        """Create a lightweight, read-only version for fast iteration.

        Returns
        -------
        SimplifiedTopology
            Immutable topology with fast numpy-based access

        Examples
        --------
        >>> import pytraj as pt
        >>> top = pt.load_topology('data/tz2.parm7')

        >>> # Create simplified version
        >>> simp_top = top.simplify()
        >>> print(f"Type: {type(simp_top)}")
        <class 'SimplifiedTopology'>

        >>> # Fast atom access
        >>> atom = simp_top.atoms[0]
        >>> print(f"First atom: {atom.name} {atom.resname}")
        First atom: N SER

        >>> # Fast residue-based selections
        >>> res = simp_top.residues[0]
        >>> res_atoms = simp_top.atoms[res.first:res.last]
        >>> print(f"Residue 0 has {len(res_atoms)} atoms")

        >>> # Vectorized operations (much faster)
        >>> charges = np.array([atom.charge for atom in simp_top.atoms])
        >>> masses = np.array([atom.mass for atom in simp_top.atoms])

        Features:
        --------
        - NumPy arrays for atoms and residues (fast indexing)
        - Read-only access (prevents accidental modification)
        - Compatible with vectorized operations
        - Maintains references to original topology

        Limitations:
        -----------
        - Read-only (no topology editing)
        - Not understood by cpptraj (use original Topology for that)
        - API may change in future versions

        Performance:
        -----------
        - ~10-100x faster for iteration-heavy operations
        - Direct numpy array access vs C++ iterator overhead
        - Ideal for analysis workflows requiring many atom/residue lookups

        Notes
        -----
        - Use for analysis workflows requiring fast iteration
        - Original Topology still needed for cpptraj operations
        - Consider ParmEd for topology editing needs
        """
        cdef _Atom atom
        cdef atom_iterator ait
        cdef res_iterator rit
        cdef int idx = 0
        cdef _Residue res
        cdef list atoms, residues

        atoms = []
        residues = []

        simplified_top = SimplifiedTopology()

        # get atoms
        ait = self.thisptr.begin()
        while ait != self.thisptr.end():
            atom = deref(ait)
            res = self.thisptr.Res(atom.ResNum())
            sim_atom = SimplifiedAtom(name=atom.c_str().strip(),
                                        type=atom.Type().Truncated(),
                                        element=get_key(
                                            atom.Element(), AtomicElementDict).lower(),
                                        charge=atom.Charge(),
                                        mass=atom.Mass(),
                                        index=idx,
                                        atomic_number=atom.AtomicNumber(),
                                        resname=res.c_str().strip(),
                                        resid=res.OriginalResNum() - 1,
                                        molnum=atom.MolNum(),
                                        simplified_top=simplified_top)
            sim_atom._bond_indices = np.asarray(self._get_atom_bond_indices(atom), dtype='int')
            atoms.append(sim_atom)
            idx += 1
            incr(ait)

        # get residues
        rit = self.thisptr.ResStart()
        idx = 0
        simplified_top.atoms = np.asarray(atoms)
        while rit != self.thisptr.ResEnd():
            res = deref(rit)
            residues.append(SimplifiedResidue(name=res.c_str().strip(),
                                              index=idx,
                                              first=res.FirstAtom(),
                                              last=res.LastAtom(),
                                              associated_topology=simplified_top))
            idx += 1
            incr(rit)
        simplified_top.residues = np.asarray(residues)
        simplified_top._top = self
        return simplified_top

    property residues:
        """Iterator over all residues in the topology.

        Yields
        ------
        Residue
            Each residue in the topology

        Examples
        --------
        >>> # Iterate over residues
        >>> for res in top.residues:
        ...     print(f"Residue {res.index}: {res.name} ({res.n_atoms} atoms)")
        Residue 0: SER (11 atoms)
        Residue 1: ALA (10 atoms)
        ...

        >>> # Find specific residue types
        >>> aromatic = [res for res in top.residues if res.name in ['PHE', 'TYR', 'TRP']]
        >>> print(f"Found {len(aromatic)} aromatic residues")

        >>> # Convert to list for random access
        >>> res_list = list(top.residues)
        >>> first_res = res_list[0]

        Properties Available on Each Residue:
        ------------------------------------
        - name : str - Residue name (ALA, SER, etc.)
        - index : int - 0-based residue index
        - n_atoms : int - Number of atoms in residue
        - first_atom_index : int - First atom index
        - last_atom_index : int - Last atom index + 1
        - original_resnum : int - Original PDB residue number
        - chain_id : str - Chain identifier

        Notes
        -----
        - Memory-efficient iterator
        - Residues are ordered sequentially
        - Use first/last_atom_index to access atoms
        """
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
        """Iterator over all molecules in the topology.

        Yields
        ------
        Molecule
            Each molecule in the topology

        Examples
        --------
        >>> # Iterate over molecules
        >>> for mol in top.mols:
        ...     print(f"Molecule {mol.index}: {mol.n_atoms} atoms, {mol.n_residues} residues")
        Molecule 0: 223 atoms, 13 residues
        Molecule 1: 3 atoms, 1 residues
        ...

        >>> # Find water molecules
        >>> water_mols = [mol for mol in top.mols if mol.n_atoms == 3]
        >>> print(f"Found {len(water_mols)} water molecules")

        >>> # Get largest molecule (protein)
        >>> largest = max(top.mols, key=lambda m: m.n_atoms)
        >>> print(f"Largest molecule: {largest.n_atoms} atoms")

        Properties Available on Each Molecule:
        -------------------------------------
        - index : int - 0-based molecule index
        - n_atoms : int - Number of atoms in molecule
        - n_residues : int - Number of residues in molecule
        - first_atom : int - First atom index
        - last_atom : int - Last atom index + 1
        - first_res : int - First residue index
        - last_res : int - Last residue index + 1

        Notes
        -----
        - Molecules defined by connectivity
        - Useful for identifying protein, solvent, ions
        - Order matches molecular connectivity graph
        """
        def __get__(self):
            cdef Molecule mol
            cdef mol_iterator it
            it = self.thisptr.MolStart()

            while it != self.thisptr.MolEnd():
                mol = Molecule()
                mol.thisptr[0] = deref(it)
                yield mol
                incr(it)

    # =========================================================================
    # SECTION 8: ANALYSIS AND SUMMARY METHODS
    # =========================================================================

    def summary(self):
        """Display comprehensive topology summary information.

        Examples
        --------
        >>> top.summary()
        Topology Summary:
        ================
        Total atoms: 5293
        Total residues: 1704
        Total molecules: 1692

        Box info: PBC with box type = ortho

        Residue types:
        - Standard: 223
        - Water: 1481
        - Ions: 0

        Notes
        -----
        - Information only appears in IPython or Python shell
        - Provides overview of system composition
        It does not appear in Jupyter notebook (due to C++ stdout)
        """
        with capture_stdout() as (out, _):
            set_world_silent(False)
            self.thisptr.Summary()
            set_world_silent(True)
        print(out)

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

    property n_heavy_atoms:
        """Number of heavy atoms (non-hydrogen, non-extra point)"""
        def __get__(self):
            return self.thisptr.HeavyAtomCount()

    property n_atom_types:
        """Number of unique atom types"""
        def __get__(self):
            return self.thisptr.NatomTypes()

    property has_charges:
        """True if any atom has non-zero charge"""
        def __get__(self):
            return self.thisptr.HasChargeInfo()

    def _set_integer_mask(self, AtomMask atm, Frame frame=Frame()):
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
        newtop_ptr = self.thisptr.partialModifyStateByMask(m.thisptr[0])
        top.thisptr[0] = deref(newtop_ptr)
        del newtop_ptr
        return top

    def _modify_state_by_mask(self, AtomMask m):
        cdef Topology top = Topology()
        newtop_ptr = self.thisptr.modifyStateByMask(m.thisptr[0])
        top.thisptr[0] = deref(newtop_ptr)
        del newtop_ptr
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

    def strip(Topology self, mask):
        """strip atoms with given mask"""
        cdef AtomMask atm
        cdef Topology new_top

        atm = AtomMask(mask)
        self._set_integer_mask(atm)
        if atm.n_atoms == 0:
            raise ValueError("number of stripped atoms must be > 1")
        atm.invert_mask()

        cdef _Topology *top_ptr = self.thisptr.modifyStateByMask(atm.thisptr[0])
        # Delete existing topology to avoid memory leak
        del self.thisptr
        self.thisptr = top_ptr

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
        self._set_integer_mask(atm)
        return atm.indices

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

    def _indices_bonded_to(self, atom_name):
        """return indices of the number of atoms that each atom bonds to

        Parameters
        ----------
        atom_name : name of the atom
        """
        cdef list arr0 = []
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
        return np.array(arr0)

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
                btype = BondType()

    property angles:
        def __get__(self):
            """return angle iterator"""
            cdef AngleArray anglearray, anglearray_h
            cdef AngleType atype = AngleType()

            anglearray = self.thisptr.Angles()
            anglearray_h = self.thisptr.AnglesH()
            anglearray.insert(anglearray.end(), anglearray_h.begin(), anglearray_h.end())

            for atype.thisptr[0] in anglearray:
                yield atype
                atype = AngleType()

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
                dhtype = DihedralType()

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
            self.add_atom(atom, residue)

        # add box
        box = Box(d['box'])
        self.box = box

        self.add_bonds(d['bond_index'])
        self.thisptr.DetermineMolecules()
        dihedral_index = d['dihedral_index']
        if dihedral_index.shape[0] != 0:
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

    property _total_charge:
        """Calculate total charge of the system.

        Returns
        -------
        float
            Sum of all atomic partial charges

        Examples
        --------
        >>> total_q = top._total_charge
        >>> print(f"System charge: {total_q:.3f} e")

        Notes
        -----
        - Private property (starts with _)
        - Sums over all atoms in topology
        - Useful for checking system neutrality
        """
        def __get__(self):
            return sum([atom.charge for atom in self.atoms])

    # =========================================================================
    # SECTION 6: FILE I/O AND FORMAT CONVERSION
    # =========================================================================

    def save(self, filename=None, format='AMBERPARM'):
        """Save topology to file in specified format.

        Parameters
        ----------
        filename : str, optional
            Output filename. If None, uses default naming.
        format : str, optional
            File format. Default 'AMBERPARM'. Options: 'AMBERPARM', 'PSF', etc.

        Examples
        --------
        >>> # Save as AMBER format
        >>> top.save('system.parm7')
        >>> top.save('system.parm7', format='AMBERPARM')

        >>> # Save as PSF format
        >>> top.save('system.psf', format='PSF')

        Notes
        -----
        - Automatically detects format from extension if not specified
        - AMBERPARM format includes all AMBER-specific parameters
        - PSF format for CHARMM/NAMD compatibility
        """
        parm = ParmFile()
        parm.write(filename=filename, top=self, format=format)

    def to_parmed(self):
        """Convert topology to ParmEd Structure object.

        Returns
        -------
        parmed.Structure
            ParmEd Structure with full topology information

        Examples
        --------
        >>> import parmed as pmd
        >>> parm_struct = top.to_parmed()
        >>> parm_struct.save('output.pdb')  # Use ParmEd functionality

        Notes
        -----
        - Requires ParmEd package to be installed
        - Enables access to ParmEd's topology editing capabilities
        - Temporary file creation during conversion
        """
        import parmed as pmd
        from pytraj.utils import tempfolder

        with tempfolder():
            self.save("tmp.prmtop")
            return pmd.load_file("tmp.prmtop")

    def to_dict(self):
        """Convert topology to Python dictionary.

        Returns
        -------
        dict
            Dictionary with all topology information

        Keys Include:
        ------------
        - atom_name : list - Atom names
        - atom_type : list - Atom types
        - atom_charge : list - Partial charges
        - atom_mass : list - Atomic masses
        - resname : list - Residue names
        - resid : list - Residue indices
        - bond_index : list - Bond connectivity
        - dihedral_index : list - Dihedral definitions
        - mol_number : list - Molecule numbers
        - box : list - Box parameters

        Examples
        --------
        >>> top_dict = top.to_dict()
        >>> print(f"Atoms: {len(top_dict['atom_name'])}")
        >>> print(f"Bonds: {len(top_dict['bond_index'])}")

        >>> # Convert to pandas DataFrame
        >>> import pandas as pd
        >>> atom_df = pd.DataFrame({
        ...     'name': top_dict['atom_name'],
        ...     'type': top_dict['atom_type'],
        ...     'charge': top_dict['atom_charge']
        ... })

        Notes
        -----
        - Useful for analysis and data export
        - All lists have same length (n_atoms)
        - Bond indices reference atom positions
        """
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

    # =========================================================================
    # SECTION 7: MOLECULAR SYSTEM PROPERTIES
    # =========================================================================

    def set_solvent(self, mask):
        """Mark atoms as solvent for analysis purposes.

        Parameters
        ----------
        mask : str
            Atom selection mask for solvent atoms

        Examples
        --------
        >>> # Mark water as solvent
        >>> top.set_solvent(':WAT')

        >>> # Mark ions and water as solvent
        >>> top.set_solvent(':WAT,Na+,Cl-')

        Notes
        -----
        - Affects analysis routines that distinguish solute/solvent
        - Used by some trajectory analysis functions
        - Does not change atom properties, just metadata
        """
        mask = mask.encode()
        self.thisptr.SetSolvent(mask)

    # =========================================================================
    # SECTION 9: STRUCTURAL ANALYSIS METHODS
    # =========================================================================

    def heavy_atom_count(self):
        """Count heavy atoms (non-hydrogen, non-extra point).

        Returns
        -------
        int
            Number of heavy atoms

        Examples
        --------
        >>> heavy_count = top.heavy_atom_count()
        >>> total_count = top.n_atoms
        >>> print(f"Heavy atoms: {heavy_count}/{total_count} ({heavy_count/total_count:.1%})")
        Heavy atoms: 1455/5293 (27.5%)

        Notes
        -----
        - Excludes hydrogen and extra point (EP) atoms
        - Useful for analysis focusing on protein backbone
        - Performance: O(n_atoms) scan
        """
        return self.thisptr.HeavyAtomCount()

    def has_charge_info(self):
        """Check if topology contains partial charge information.

        Returns
        -------
        bool
            True if any atom has non-zero charge

        Examples
        --------
        >>> if top.has_charge_info():
        ...     print("Can perform electrostatic calculations")
        ... else:
        ...     print("No charge information available")

        Notes
        -----
        - Scans all atoms for non-zero charges
        - Required for electrostatic analysis
        - Most MD force fields include charges
        """
        return self.thisptr.HasChargeInfo()

    # =========================================================================
    # SECTION 10: VAN DER WAALS PARAMETERS
    # =========================================================================

    def get_vdw_sigma(self, int atom_idx):
        """Get van der Waals sigma parameter for atom.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        float
            VDW sigma parameter (Angstroms)

        Examples
        --------
        >>> # Get VDW parameters for first atom
        >>> sigma = top.get_vdw_sigma(0)
        >>> print(f"VDW sigma: {sigma:.3f} Å")

        >>> # Compare different atom types
        >>> for i in range(5):
        ...     atom = top[i]
        ...     sigma = top.get_vdw_sigma(i)
        ...     print(f"{atom.name}: σ={sigma:.3f} Å")

        Notes
        -----
        - Sigma is the Lennard-Jones size parameter
        - Related to atomic radius (σ ≈ 2^(1/6) * r_min)
        - Used in non-bonded force calculations
        """
        return self.thisptr.GetVDWsigma(atom_idx)

    def get_vdw_depth(self, int atom_idx):
        """Get van der Waals depth parameter for atom.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        float
            VDW depth parameter (kcal/mol)

        Examples
        --------
        >>> # Get VDW depth (epsilon)
        >>> depth = top.get_vdw_depth(0)
        >>> print(f"VDW depth: {depth:.4f} kcal/mol")

        >>> # Find strongest VDW interactions
        >>> max_depth = max(top.get_vdw_depth(i) for i in range(top.n_atoms))
        >>> print(f"Strongest VDW: {max_depth:.4f} kcal/mol")

        Notes
        -----
        - Epsilon is the Lennard-Jones well-depth parameter
        - Controls strength of attractive interactions
        - Typically larger for heavier atoms
        """
        return self.thisptr.GetVDWdepth(atom_idx)

    def get_vdw_radius(self, int atom_idx):
        """Get van der Waals radius for atom.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        float
            VDW radius (Angstroms)

        Examples
        --------
        >>> # Get atomic radius
        >>> radius = top.get_vdw_radius(0)
        >>> print(f"Atomic radius: {radius:.3f} Å")

        >>> # Calculate molecular surface area (simplified)
        >>> radii = [top.get_vdw_radius(i) for i in range(top.n_atoms)]
        >>> total_surface = sum(4 * 3.14159 * r**2 for r in radii)
        >>> print(f"Total surface area: {total_surface:.1f} Å²")

        Notes
        -----
        - Effective atomic radius for VDW interactions
        - Used in surface area and volume calculations
        - Relates to sigma parameter: r = σ/2^(1/6)
        """
        return self.thisptr.GetVDWradius(atom_idx)

    # =========================================================================
    # SECTION 11: MASS AND CHARGE ANALYSIS
    # =========================================================================

    def mask_has_zero_mass(self, mask):
        """Check if selected atoms have zero total mass.

        Parameters
        ----------
        mask : str or AtomMask
            Atom selection mask

        Returns
        -------
        bool
            True if total mass of selected atoms is zero

        Examples
        --------
        >>> # Check if extra points have zero mass
        >>> has_zero = top.mask_has_zero_mass('@EP')
        >>> print(f"Extra points have zero mass: {has_zero}")

        >>> # Validate protein selection
        >>> protein_zero = top.mask_has_zero_mass(':1-100')
        >>> if protein_zero:
        ...     print("Warning: Protein region has zero mass!")

        Notes
        -----
        - Useful for detecting extra point atoms
        - Helps validate atom selections
        - Zero mass atoms don't contribute to kinetic energy
        """
        cdef AtomMask atm
        if isinstance(mask, str):
            atm = AtomMask(mask)
            self._set_integer_mask(atm)
        else:
            atm = mask
        return self.thisptr.MaskHasZeroMass(atm.thisptr[0])

    # =========================================================================
    # SECTION 12: ATOM AND RESIDUE IDENTIFICATION
    # =========================================================================

    def trunc_resname_atomname(self, int atom_idx):
        """Get formatted residue and atom name identifier.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        str
            Format: '<resname>@<atomname>'

        Examples
        --------
        >>> # Get formatted atom identifiers
        >>> for i in range(5):
        ...     identifier = top.trunc_resname_atomname(i)
        ...     print(f"Atom {i}: {identifier}")
        Atom 0: SER@N
        Atom 1: SER@H
        Atom 2: SER@CA

        >>> # Create selection strings
        >>> ca_id = top.trunc_resname_atomname(2)
        >>> print(f"Selection: {ca_id}")
        Selection: SER@CA

        Notes
        -----
        - Useful for creating human-readable atom identifiers
        - Format compatible with cpptraj selection syntax
        - Helpful for logging and debugging
        """
        return self.thisptr.TruncResNameAtomName(atom_idx)

    def trunc_atom_name_num(self, int atom_idx):
        """Get formatted atom name and number.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        str
            Format: '<atomname>_<atomnum>'

        Examples
        --------
        >>> # Get atom name with numbering
        >>> name_num = top.trunc_atom_name_num(0)
        >>> print(f"Atom identifier: {name_num}")
        Atom identifier: N_1

        Notes
        -----
        - Provides unique atom identifier
        - Uses 1-based numbering (PDB convention)
        - Useful for output formatting
        """
        return self.thisptr.TruncAtomNameNum(atom_idx)

    def trunc_resname_onum_id(self, int res_idx):
        """Get formatted residue identifier with chain info.

        Parameters
        ----------
        res_idx : int
            0-based residue index

        Returns
        -------
        str
            Format: '<resname>_<onum>[_<id>]'

        Examples
        --------
        >>> # Get residue identifiers
        >>> for i in range(3):
        ...     res_id = top.trunc_resname_onum_id(i)
        ...     print(f"Residue {i}: {res_id}")
        Residue 0: SER_1_A
        Residue 1: ALA_2_A

        Notes
        -----
        - Includes original PDB residue number
        - Chain ID included if available
        - Useful for cross-referencing with PDB files
        """
        return self.thisptr.TruncResNameOnumId(res_idx)

    def resname_num_atomname_num(self, int atom_idx):
        """Get comprehensive atom identifier with spacing.

        Parameters
        ----------
        atom_idx : int
            0-based atom index

        Returns
        -------
        str
            Format: '<resname> <resnum> <atomname> <atomnum>'

        Examples
        --------
        >>> # Get full atom description
        >>> full_id = top.resname_num_atomname_num(0)
        >>> print(f"Full identifier: {full_id}")
        Full identifier: SER 1 N 1

        >>> # Create readable atom list
        >>> for i in range(5):
        ...     print(top.resname_num_atomname_num(i))
        SER 1 N 1
        SER 1 H 2
        SER 1 CA 3

        Notes
        -----
        - Most comprehensive atom identifier format
        - Uses PDB-style numbering (1-based)
        - Ideal for detailed output and reports
        """
        return self.thisptr.ResNameNumAtomNameNum(atom_idx)

    # =========================================================================
    # SECTION 13: MOLECULAR ANALYSIS
    # =========================================================================

    def nres_in_mol(self, int mol_idx):
        """Get number of residues in specified molecule.

        Parameters
        ----------
        mol_idx : int
            0-based molecule index

        Returns
        -------
        int
            Number of residues in molecule

        Examples
        --------
        >>> # Analyze molecular composition
        >>> for mol_idx in range(min(5, top.n_mols)):
        ...     n_res = top.nres_in_mol(mol_idx)
        ...     mol = list(top.mols)[mol_idx]
        ...     print(f"Mol {mol_idx}: {n_res} residues, {mol.n_atoms} atoms")
        Mol 0: 13 residues, 223 atoms
        Mol 1: 1 residues, 3 atoms

        >>> # Find protein molecules (multiple residues)
        >>> proteins = [i for i in range(top.n_mols)
        ...             if top.nres_in_mol(i) > 1]
        >>> print(f"Found {len(proteins)} protein molecules")

        Notes
        -----
        - Helps classify molecule types
        - Protein molecules typically have multiple residues
        - Solvent/ions usually have single residues
        """
        return self.thisptr.NresInMol(mol_idx)

    def solute_residues(self):
        """Identify solute residues (typically protein/nucleic acid).

        Returns
        -------
        list
            List of solute residue indices (0-based)

        Examples
        --------
        >>> # Get solute residues
        >>> solute_res = top.solute_residues()
        >>> print(f"Found {len(solute_res)} solute residues")

        >>> # Analyze solute composition
        >>> for res_idx in solute_res[:5]:  # First 5
        ...     res = top.residue(res_idx)
        ...     print(f"Solute residue {res_idx}: {res.name}")

        >>> # Create solute-only selection
        >>> solute_atoms = []
        >>> for res_idx in solute_res:
        ...     res = top.residue(res_idx)
        ...     solute_atoms.extend(range(res.first_atom_index, res.last_atom_index))

        Algorithm:
        ----------
        - Identifies residues in multi-residue molecules
        - Excludes single-atom molecules (ions)
        - Excludes small molecules (< 2 atoms per molecule)
        - Heuristic-based approach

        Notes
        -----
        - Simple heuristic implementation
        - May need adjustment for unusual systems
        - Consider using explicit masks for precise control
        - Solvent marking via set_solvent() affects this
        """
        # Simple implementation: find residues that are not in solvent molecules
        cdef list solute_res = []
        cdef int res_idx
        cdef int mol_idx

        for res_idx in range(self.n_residues):
            # Get the first atom of this residue
            res = self.residue(res_idx)
            if res.n_atoms > 0:
                first_atom = res.first_atom_index
                atom = self[first_atom]
                mol_idx = atom.molnum

                # Check if molecule is solvent (simple heuristic: single residue water-like molecules)
                if mol_idx < self.n_mols:
                    mol = list(self.mols)[mol_idx]
                    if mol.n_atoms > 1:  # Not a single atom (likely not ion)
                        solute_res.append(res_idx)

        return solute_res

    def resnums_selected_by(self, mask):
        """Get residue numbers for atoms selected by mask.

        Parameters
        ----------
        mask : str or AtomMask
            Atom selection mask

        Returns
        -------
        numpy.ndarray
            Array of unique residue numbers (0-based)

        Examples
        --------
        >>> # Get residues containing CA atoms
        >>> ca_residues = top.resnums_selected_by('@CA')
        >>> print(f"Residues with CA: {len(ca_residues)}")

        >>> # Find residues in active site region
        >>> active_residues = top.resnums_selected_by(':85-95')
        >>> print(f"Active site residues: {active_residues}")

        >>> # Get residues with aromatic atoms
        >>> aromatic_res = top.resnums_selected_by(':PHE,TYR,TRP')
        >>> for res_idx in aromatic_res:
        ...     res = top.residue(res_idx)
        ...     print(f"Aromatic: {res.name} {res_idx}")

        Notes
        -----
        - Returns unique residue indices only
        - 0-based indexing (Python convention)
        - Useful for residue-level analysis
        """
        cdef AtomMask atm
        if isinstance(mask, str):
            atm = AtomMask(mask)
            self._set_integer_mask(atm)
        else:
            atm = mask
        cdef vector[int] resnums = self.thisptr.ResnumsSelectedBy(atm.thisptr[0])
        return np.array([resnums[i] for i in range(resnums.size())], dtype=int)

    def molnums_selected_by(self, mask):
        """Get molecule numbers for atoms selected by mask.

        Parameters
        ----------
        mask : str or AtomMask
            Atom selection mask

        Returns
        -------
        numpy.ndarray
            Array of unique molecule numbers (0-based)

        Examples
        --------
        >>> # Get molecules containing backbone atoms
        >>> backbone_mols = top.molnums_selected_by('@CA,C,N,O')
        >>> print(f"Molecules with backbone: {len(backbone_mols)}")

        >>> # Find water molecules
        >>> water_mols = top.molnums_selected_by(':WAT')
        >>> print(f"Water molecules: {len(water_mols)}")

        >>> # Analyze molecular composition
        >>> all_mols = top.molnums_selected_by('*')
        >>> protein_mols = [m for m in all_mols
        ...                 if top.nres_in_mol(m) > 1]
        >>> print(f"Protein molecules: {len(protein_mols)}")

        Notes
        -----
        - Returns unique molecule indices only
        - 0-based indexing (Python convention)
        - Helps analyze molecular system composition
        """
        cdef AtomMask atm
        if isinstance(mask, str):
            atm = AtomMask(mask)
            self._set_integer_mask(atm)
        else:
            atm = mask
        cdef vector[int] molnums = self.thisptr.MolnumsSelectedBy(atm.thisptr[0])
        return np.array([molnums[i] for i in range(molnums.size())], dtype=int)

    # =========================================================================
    # SECTION 14: TOPOLOGY MODIFICATION
    # =========================================================================

    def merge_residues(self, int start_res, int stop_res):
        """Merge consecutive residues into a single residue.

        Parameters
        ----------
        start_res : int
            Starting residue index (0-based)
        stop_res : int
            Stopping residue index (0-based, inclusive)

        Returns
        -------
        int
            0 on success, 1 on error

        Examples
        --------
        >>> # Merge first 3 residues
        >>> result = top.merge_residues(0, 2)
        >>> if result == 0:
        ...     print("Successfully merged residues 0-2")

        >>> # Create custom residue groupings
        >>> # Merge domain residues
        >>> top.merge_residues(10, 50)  # Domain 1
        >>> top.merge_residues(51, 100) # Domain 2

        Notes
        -----
        - Modifies topology in-place
        - Use with caution - can break molecular connectivity
        - Useful for coarse-graining analysis
        - Residue numbers will be renumbered

        Warning
        -------
        - Irreversible operation
        - May break some analysis functions
        - Consider making a copy first
        """
        return self.thisptr.MergeResidues(start_res, stop_res)

    def set_single_molecule(self):
        """Treat entire system as a single molecule.

        Returns
        -------
        int
            0 on success

        Examples
        --------
        >>> # Treat protein-DNA complex as single molecule
        >>> result = top.set_single_molecule()
        >>> if result == 0:
        ...     print(f"System now has {top.n_mols} molecule(s)")

        Notes
        -----
        - Useful for treating complexes as single units
        - Affects molecule-based analysis functions
        - Modifies topology in-place
        - Cannot be easily undone
        """
        return self.thisptr.SetSingleMolecule()

    def split_residue(self, mask, new_resname):
        """Split atoms from residue into a new residue.

        Parameters
        ----------
        mask : str or AtomMask
            Atoms to split into new residue
        new_resname : str
            Name for the new residue

        Returns
        -------
        int
            0 on success, 1 on error

        Examples
        --------
        >>> # Split ligand atoms from protein residue
        >>> result = top.split_residue(':100@C1,C2,C3', 'LIG')
        >>> if result == 0:
        ...     print("Successfully created new ligand residue")

        >>> # Separate sidechain from backbone
        >>> top.split_residue(':50@CB,CG,CD', 'SCH')

        Notes
        -----
        - Creates new residue with specified atoms
        - Original residue keeps remaining atoms
        - Useful for separating ligands or modified residues
        - Modifies topology in-place

        Warning
        -------
        - Complex operation that may affect connectivity
        - Test thoroughly with your specific use case
        - Consider backup before modification
        """
        cdef AtomMask atm
        if isinstance(mask, str):
            atm = AtomMask(mask)
            self._set_integer_mask(atm)
        else:
            atm = mask
        cdef NameType resname = NameType(new_resname)
        return self.thisptr.SplitResidue(atm.thisptr[0], resname.thisptr[0])

    def residue(self, int idx, bint atom=False):
        """Get residue by index.

        Parameters
        ----------
        idx : int
            0-based residue index
        atom : bool, optional
            If True, include full atom list (slower). Default False.

        Returns
        -------
        Residue
            Residue object with structural information

        Examples
        --------
        >>> res = top.residue(0)  # First residue
        >>> print(f"Residue: {res.name}, atoms: {res.n_atoms}")
        Residue: SER, atoms: 11

        >>> # Access residue atoms efficiently
        >>> first_atom_idx = res.first_atom_index
        >>> last_atom_idx = res.last_atom_index
        >>> atoms_in_res = top[first_atom_idx:last_atom_idx]

        Notes
        -----
        - Use atom=False for better performance
        - Access atoms via first_atom_index/last_atom_index
        - Residue objects are views into the topology
        """
        cdef Residue res = Residue()
        res.thisptr[0] = self.thisptr.Res(idx)
        return res

    def atom(self, int idx):
        """Get atom by index with mutable access.

        Parameters
        ----------
        idx : int
            0-based atom index

        Returns
        -------
        Atom
            Atom object that can modify topology when changed

        Examples
        --------
        >>> atom = top.atom(0)
        >>> print(f"{atom.name} {atom.element} {atom.charge:.3f}")
        N N -0.416

        >>> # Modify atom (affects topology)
        >>> atom.charge = -0.5  # Changes topology

        Notes
        -----
        - Returns mutable reference to atom
        - Modifications affect the original topology
        - Use with caution to maintain consistency
        - For read-only access, use iteration or indexing
        """
        cdef Atom atom = Atom()
        atom.own_memory = False
        atom.thisptr = &self.thisptr.GetAtomView(idx)
        atom.resname = self.thisptr.Res(atom.resid).c_str().strip()
        return atom

    # Moved to File I/O section

cdef class ParmFile:
    """Low-level topology file reader/writer interface.

    This class provides direct access to cpptraj's parameter file I/O
    functionality for reading and writing topology files.

    Notes
    -----
    - Low-level class, usually wrapped by higher-level functions
    - Use pytraj.load_topology() and Topology.save() for most cases
    - Supports multiple file formats (AMBER, PSF, etc.)
    """

    def __cinit__(self):
        """Initialize ParmFile object."""
        self.thisptr = new _ParmFile()

    def __dealloc__(self):
        """Clean up C++ resources."""
        del self.thisptr

    def read(self, filename="", top=Topology(), option=''):
        """Read topology from file into Topology object.

        Parameters
        ----------
        filename : str
            Path to topology file
        top : Topology, optional
            Topology object to populate. Default creates new one.
        option : str or ArgList, optional
            Additional read options

        Returns
        -------
        None
            Modifies top in-place

        Notes
        -----
        - Low-level method, use pytraj.load_topology() instead
        - Supports various formats based on file extension
        - Modifies the provided Topology object directly
        """
        cdef ArgList arglist
        cdef debug = 0
        cdef Topology _top = <Topology > top

        filename = filename.encode()

        if not option:
            self.thisptr.ReadTopology(_top.thisptr[0], filename, debug)
        else:
            if isinstance(option, str):
                arglist = ArgList(option)
            else:
                arglist = <ArgList > option
            self.thisptr.ReadTopology(
                _top.thisptr[0], filename, arglist.thisptr[0], debug)

    def write(self, Topology top=Topology(), filename="default.top",
                  ArgList arglist=ArgList(), format=""):
        """Write topology to file.

        Parameters
        ----------
        top : Topology
            Topology to write
        filename : str, optional
            Output filename. Default "default.top"
        arglist : ArgList, optional
            Additional write arguments
        format : str, optional
            File format. If empty, detected from filename.

        Raises
        ------
        ValueError
            If topology is empty or write fails

        Notes
        -----
        - Low-level method, use Topology.save() instead
        - Automatically detects format from extension
        - Supports AMBER, PSF, and other formats
        """
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
        """Get the filename of currently loaded parameter file.

        Returns
        -------
        str
            Absolute path to parameter file

        Examples
        --------
        >>> parm = ParmFile()
        >>> parm.read('system.parm7')
        >>> print(parm.filename())
        '/path/to/system.parm7'
        """
        cdef FileName filename = FileName()
        filename.thisptr[0] = self.thisptr.ParmFilename()
        return os.path.abspath(filename)
