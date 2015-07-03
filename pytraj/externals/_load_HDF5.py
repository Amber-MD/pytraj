from __future__ import absolute_import
from .six import string_types


def load_hdf5(filename_or_buffer, autoconvert=True, restype=None, top=None):
    """"load hd5f format from openmm (?)

    Parameters
    ---------
    filename_or_buffer : str or buffer
    autoconvert : bool
        if 'True' (default): convert from `nm` to `angstrom`
        if 'False': No convert

    Examples
    -------
        from mdtraj.testing import get_fn
        fname = get_fn("frame0.h5")
        import pytraj.io as io
        traj = io.load_hdf5(fname, autoconvert=False)
        print (traj)
    """

    try:
        import h5py
    except ImportError:
        raise ImportError("require h5py, HDF5 lib and numpy")

    if isinstance(filename_or_buffer, string_types):
        fh = h5py.File(filename_or_buffer, 'r')
        should_be_closed = True
    else:
        fh = filename_or_buffer
        should_be_closed = False

    traj = _load_hdf5_from_buffer(
        fh, autoconvert=autoconvert, restype=restype, top=top)

    if should_be_closed:
        fh.close()

    return traj


def _load_hdf5_from_buffer(fh, autoconvert=True, restype=None, top=None):
    import json
    from ..Topology import Topology
    from ..Trajectory import Trajectory
    from ..core import Atom, Box
    from ..core import mass_atomic_number_dict, mass_element_dict
    # NOTE: always use `np.float64` in pytraj
    if autoconvert:
        UNIT = 10.
    else:
        UNIT = 1.

    try:
        import h5py
        import numpy as np
    except ImportError:
        raise ImportError("require h5py, HDF5 lib and numpy")

    try:
        cell_lengths = fh['cell_lengths'].value * UNIT
        box_arr = np.hstack(
            (cell_lengths, fh['cell_angles'])).astype(np.float64)
        has_box = True
    except:
        has_box = False

    crd = fh['coordinates'].value.astype(np.float64)
    n_frames, n_atoms, _ = crd.shape
    if autoconvert:
        crd = crd * UNIT

    if restype is None:
        farray = Trajectory()
        farray._allocate(n_frames, n_atoms)
    elif restype == 'api.Trajectory':
        from pytraj import api
        farray = api.Trajectory()

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
                    # no infor about atom type in .h5 file from openmm (?)
                    atype = aname
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
                                mass = mass_atomic_number_dict[
                                    atom['atomic_number']]
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
    farray.top = _top

    # update coords
    if restype is None:
        farray.update_xyz(crd)
        if has_box:
            for idx, arr in enumerate(crd):
                farray[idx].box = box_arr[idx]  # auto-cast
    else:
        farray.xyz = crd

    return farray
