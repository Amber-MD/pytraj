"""Load external parm object
"""
from __future__ import absolute_import
from pytraj.utils import has_, _import_numpy
from pytraj.Trajectory import Trajectory
from pytraj.Topology import Topology
from pytraj.core.Atom import Atom
from pytraj.Frame import Frame
from pytraj.utils.check_and_assert import is_mdtraj, is_mdanalysis
_, np = _import_numpy()

__all__ = ['load_pseudo_parm']

def load_pseudo_parm(parm, guess_bond=True):
    """load_external's parm objects

    Parameters
    ---------
    parm : external Topology/Parm objects (mdtraj, parmed) 
        or Universe object (MDAnalysis)
    """
    from pytraj.core import Box
    farray = Trajectory()

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
            # NOTE: need to convert to string for some Topology (PSF, ...)
            res = atom.residue
            aname = str(atom.name)
            resname = str(res.name)

            if is_mdtraj(parm):
                atype = str(atom.name) # mdtraj
                resid = res.index
                mass = atom.element.mass
                #charge = atom.element.charge
                charge = 0.0
            elif is_mdanalysis(parm):
                # in MDAnalysis, atom.type is `int`
                atype = str(atom.type) 
                resid = atom.resid
                charge = atom.charge
                mass = atom.mass
            else:
                atype = str(atom.type) # parmed
                resid = res.idx
                charge = atom.charge
                mass = atom.mass
            atom = Atom(aname, atype, charge, mass)
            pseudotop.add_atom(atom=atom, resid=resid, resname=resname)

    if is_mdtraj(parm):
        # not sure how to get angles, dihedrals quickly
        try:
            pseudotop.add_bonds(np.array([(a.index, b.index) for (a, b) in parm.bonds]))
        except:
            pass
        # load Box in _load_mdtraj since Box is stored in traj
    elif is_mdanalysis(parm):
        # turn-off. need to check MDAnalysis
        # ack, lots of try and except
        try:
            pseudotop.add_bonds(np.asarray(parm.universe.bonds.to_indices()))
        except:
            pass

        try:
            pseudotop.add_angles(np.asarray(parm.universe.angles.to_indices()))
        except:
            pass

        try:
            pseudotop.add_dihedrals(np.asarray(parm.universe.torsions.to_indices()))
        except:
            pass

        try:
            pseudotop.box = Box(parm.dimensions.astype(np.float64))
        except:
            pass
    else:
        # TODO : add bonds, dihedrals, angles for ParmEd
        # parmed
        # add dihedrals
        try:
            bond_list = [(x.atom1.idx, x.atom2.idx)
                        for x in parm.bonds]
            angle_list = [(x.atom1.idx, x.atom2.idx, x.atom3.idx)
                        for x in parm.angles]
            dihedral_list= [(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx)
                      for x in parm.dihedrals]
            pseudotop.add_bonds(np.asarray(bond_list))
            pseudotop.add_angles(np.asarray(angle_list))
            pseudotop.add_dihedrals(np.asarray(dihedral_list))
        except:
            pass

        try:
            pseudotop.box = Box(np.array(parm.box))
        except:
            # no box
            pseudotop.box = Box()
    try:
        pseudotop._bonds_ndarray[0]
    except IndexError:
        # TODO, FIXME: really slow if loading parmed
        # >>> import parmed as pmd
        # >>> pdb = pmd.download_PDB("1o15")
        # >>> top = io.load_pseudo_parm(pdb) $ 46 s. Really?
        if guess_bond:
            pseudotop.guess_bond()
        pass
    return pseudotop
