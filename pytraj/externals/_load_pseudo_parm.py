"""Load external parm object

"""
from __future__ import absolute_import
from pytraj.utils import has_, _import_numpy
from pytraj.FrameArray import FrameArray
from pytraj.Topology import Topology
from pytraj.core.Atom import Atom
from pytraj.Frame import Frame
from pytraj.utils.check_and_assert import is_mdtraj, is_word_in_class_name

# not sure if we need this `load_mdtraj` since cpptraj can do anything :D
# might need to move to Cython level for faster loading

_, np = _import_numpy()

def load_pseudo_parm(parm):
    # TODO: fill me
    """load_external's parm objects

    Parameters
    ---------
    parm : external Topology/Parm objects (mdtraj, chemistry) 
        or Universe object (MDAnalysis)t
    """
    farray = FrameArray()

    # convert to pseudo-topology
    # to fully use Topology object in pytraj, we can do:
    # >>> farray.top = Topology(top_name) # or
    # >>> from pytraj import io
    # >>> farray.top = io.load(top_name) 

    pseudotop = Topology()
    for atom in parm.atoms:
        res = atom.residue
        aname = atom.name
        resname = res.name

        if is_mdtraj(parm):
            atype = atom.name # mdtraj
            resid = res.index
        elif is_word_in_class_name(parm, 'Universe'):
            # in MDAnalysis, atom.type is `int`
            atype = str(atom.type) 
            resid = atom.resid
        else:
            atype = atom.type # parmed
            resid = res.idx
        atom = Atom(aname, atype)
        pseudotop.add_atom(atom=atom, resid=resid, resname=resname)
    if is_word_in_class_name(parm, 'Universe'):
        pseudotop.add_bonds(np.asarray(its_obj.bonds.to_indices()))
        pseudotop.add_angles(np.asarray(its_obj.angles.to_indices()))
        pseudotop.add_dihedrals(np.asarray(its_obj.torsions.to_indices()))
    elif is_mdtraj():
        # not sure how to get bonds, angles, dihedrals quickly
        pass
    # if is ParmEd: use load_full_ParmEd()
    return pseudotop
