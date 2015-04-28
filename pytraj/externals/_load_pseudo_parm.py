"""Load external parm object
"""
from __future__ import absolute_import
from pytraj.utils import has_, _import_numpy
from pytraj.FrameArray import FrameArray
from pytraj.Topology import Topology
from pytraj.core.Atom import Atom
from pytraj.Frame import Frame
from pytraj.utils.check_and_assert import is_mdtraj, is_mdanalysis
_, np = _import_numpy()

__all__ = ['load_pseudo_parm']

def load_pseudo_parm(parm):
    """load_external's parm objects

    Parameters
    ---------
    parm : external Topology/Parm objects (mdtraj, chemistry) 
        or Universe object (MDAnalysis)t
    """
    from pytraj.core import Box
    farray = FrameArray()

    # convert to pseudo-topology
    # to fully use Topology object in pytraj, we can do:
    # >>> farray.top = Topology(top_name) # or
    # >>> from pytraj import io
    # >>> farray.top = io.load(top_name) 
    if is_mdanalysis(parm):
        #chains = parm.fragments
        chains = parm.segments
    elif is_mdtraj(parm):
        chains = parm.chains
    else:
        chains = [parm,] # fake

    pseudotop = Topology()
    for chain in chains:
        pseudotop.start_new_mol()
        for atom in chain.atoms:
            res = atom.residue
            aname = atom.name
            resname = res.name

            if is_mdtraj(parm):
                atype = atom.name # mdtraj
                resid = res.index
            elif is_mdanalysis(parm):
                # in MDAnalysis, atom.type is `int`
                atype = str(atom.type) 
                resid = atom.resid
            else:
                atype = atom.type # parmed
                resid = res.idx
            atom = Atom(aname, atype)
            # TODO : add mass too
            pseudotop.add_atom(atom=atom, resid=resid, resname=resname)

    if is_mdanalysis(parm):
        pseudotop.add_bonds(np.asarray(parm.bonds.to_indices()))
        pseudotop.add_angles(np.asarray(parm.angles.to_indices()))
        pseudotop.add_dihedrals(np.asarray(parm.torsions.to_indices()))
        pseudotop.box = Box(parm.dimensions.astype(np.float64))
    elif is_mdtraj(parm):
        # not sure how to get angles, dihedrals quickly
        pseudotop.add_bonds(np.array([(a.index, b.index) for (a, b) in parm.bonds]))
        # load Box in _load_mdtraj since Box is stored in traj
    else:
        # parmed
        pseudotop.box = Box(np.array(parm.box))
    return pseudotop
