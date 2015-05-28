from __future__ import absolute_import
from pytraj.externals.six import add_metaclass

from .. utils.check_and_assert import is_int
from .. Trajectory import Trajectory
from .. Frame import Frame
from .. core import Atom, Box
from .. Topology import Topology
from .. externals import load_pseudo_parm
from .. _shared_methods import my_str_method
from .. Topology import Topology
from . TrajectoryBaseIterator import TrajectoryBaseIterator

class TrajectoryHDF5Iterator(TrajectoryBaseIterator):
    def __init__(self, filename, top=None):
        self._top = Topology()
        self._filename = filename
        self._fh = None
        self._box_arr = None

        if self._filename is not None:
            self.load(filename, top)

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        frame = Frame(self._n_atoms)
        for idx, xyz in enumerate(self._crd):
            frame._fast_copy_from_xyz(xyz)
            if self.has_box():
                frame.box = Box(self._box_arr[idx])
            yield frame

    def __getitem__(self, idx):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._fh.close()

    @property
    def n_frames(self):
        return self._n_frames

    @property
    def size(self):
        return self.n_frames

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def xyz(self):
        return self._crd

    @property
    def filename(self):
        return self._filename

    @property
    def top(self):
        return self._top

    @top.setter
    def top(self, value):
        self._top = value

    def has_box(self):
        return self._has_box

    def load(self, filename, top=None, mode='r'):
        try:
            import h5py
            import numpy as np
            import json
        except ImportError:
            raise ImportError("requilre h5py to read HDF5 file")
        fh = h5py.File(filename, mode)
        self._fh = fh

        try:
            cell_lengths = fh['cell_lengths'].value
            box_arr = np.hstack((cell_lengths, fh['cell_angles'])).astype(np.float64)
            has_box = True
            self._box_arr = box_arr
        except:
            has_box = False

        crd = fh['coordinates'].value.astype('f8')
        n_frames, n_atoms, _ = crd.shape
        self._n_frames = n_frames
        self._n_atoms = n_atoms
        self._crd = crd
        self._filename = fh.filename

        # create Topology
        if top is not None:
            _top = top
        else:
            top_txt = fh['topology']
            h5_topology = json.loads(top_txt.value.tostring().decode())
            _top = Topology()
            for chain in h5_topology['chains']:
                _top.start_new_mol()
                for residue in chain['residues']:
                    resname = residue['name']
                    resid = residue['index']
                    for atom in residue['atoms']:
                        aname = atom['name']
                        atype = aname # no infor about atom type in .h5 file from openmm (?)
                        try:
                            charge = atom['charge']
                        except:
                            charge = 0.0
                        try:
                            mass = atom['mass']
                        except:
                            try:
                                mass = mass_element_dict[atom['element']]
                            except:
                                try:
                                    mass = mass_atomic_number_dict[atom['atomic_number']]
                                except:
                                    mass = 1.0
                        atom = Atom(aname, atype, charge, mass)
                        _top.add_atom(atom=atom, resid=resid, resname=resname)
            # add bonds
            # Note: no PBC info for top
            _top.add_bonds(np.asarray(h5_topology['bonds']))
            # naively assigne box info from 1st frame
            if has_box:
                _top.box = Box(box_arr[0])
                self._has_box = True
            else:
                self._has_box= False

        self.top = _top
