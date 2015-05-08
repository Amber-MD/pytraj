from __future__ import absolute_import
from ..Topology import Topology
from ..Trajectory import Trajectory
from ..Frame import Frame
from ..core import Atom, Box

def load_hd5f(filename, autoconvert=True):
    """"load hd5f format from openmm (?)

    Parameters
    ---------
    filename : str 
    autoconvert : bool
        if 'True' (default): convert from `nm` to `angstrom`
        if 'False': No convert

    Examples
    -------
        from mdtraj.testing import get_fn
        fname = get_fn("frame0.h5")
        import pytraj.io as io
        traj = io.load_hd5f(fname, autoconvert=False)
        print (traj)
    """
    import json
    # NOTE: always use `np.float64` in pytraj
    if autoconvert:
        UNIT = 10.
    else:
        UNIT = 1.

    try:
        import h5py
        import numpy as np
    except ImportError:
        raise ImportError("require h5py, HD5F lib and numpy")

    fh = h5py.File(filename, 'r')
    crd = fh['coordinates']
    try:
        cell_lengths = fh['cell_lengths'].value * UNIT
        box_arr = np.hstack((cell_lengths, fh['cell_angles'])).astype(np.float64)
        has_box = True
    except:
        has_box = False
    shape = crd.shape

    farray = Trajectory()
    # create Topology
    top_txt = fh['topology']
    h5_topology = json.loads(top_txt.value.tostring().decode())
    top = Topology()
    for chain in h5_topology['chains']:
        top.start_new_mol()
        for residue in chain['residues']:
            resname = residue['name']
            resid = residue['index']
            for atom in residue['atoms']:
                aname = atom['name']
                atype = aname # no infor about atom type in .h5 file from openmm (?)
                atom = Atom(aname, atype)
                top.add_atom(atom=atom, resid=resid, resname=resname)
    # add bonds
    # Note: no PBC info for top
    top.add_bonds(np.asarray(h5_topology['bonds']))
    # naively assigne box info from 1st frame
    if has_box:
        top.box = Box(box_arr[0])
    farray.top = top

    # update coords
    farray.resize(shape[0])
    for idx, arr in enumerate(crd):
        # allocate memory for Frame
        farray[idx] = Frame(shape[1])
        # update coords by cython memoryview method (faster)
        farray[idx, :] = arr.astype(np.float64) * UNIT
        if has_box:
            farray[idx].box = Box(box_arr[idx])

    return farray
