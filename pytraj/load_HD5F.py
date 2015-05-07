from __future__ import absolute_import
from .Topology import Topology
from .Trajectory import Trajectory
from .Frame import Frame

def load_hd5f(filename):
    farray = Trajectory()
    try:
        import h5py
        import numpy as np
    except ImportError:
        raise ImportError("require h5py, HD5F lib and numpy")

    fh = h5py.File(filename, 'r')
    crd = fh['coordinates']
    shape = crd.shape

    farray.resize(shape[0])
    for idx, arr in enumerate(crd):
        farray[idx] = Frame(shape[1])
        farray[idx].set_from_crd(arr.flatten().astype(np.float64))

    # TODO : load Topology from h5py object
    #pseudotop = Topology()
    #for mdatom in mtop.atoms:
    #    atom = Atom(mdatom.name, mdatom.name)
    #    mdres = mdatom.residue
    #    pseudotop.add_atom(atom=atom, resid=mdres.index, resname=mdres.name)
    #return farray
