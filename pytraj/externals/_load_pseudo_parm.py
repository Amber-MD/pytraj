"""Load external parm object
"""
from __future__ import absolute_import

__all__ = ['load_pseudo_parm']


def load_pseudo_parm(parm, guess_bond=True):
    """load_external's parm objects

    Parameters
    ---------
    parm : external Topology/Parm objects (mdtraj, parmed) 
        or Universe object (MDAnalysis)
    """
    from pytraj.utils import _import_numpy
    from pytraj.Trajectory import Trajectory
    from pytraj.Topology import Topology
    from pytraj.core.Atom import Atom
    from pytraj.utils.check_and_assert import is_mdtraj, is_mdanalysis
    _, np = _import_numpy()

    from pytraj.core import Box

    if is_mdanalysis(parm):
        chains = parm.segments
    elif is_mdtraj(parm):
        chains = parm.chains
    else:
        chains = [parm,] # fake

    pseudotop = Topology()
    i = 0
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
                if atom.type:
                    atype = str(atom.type) # parmed
                else:
                    atype = atom.name
                resid = res.idx
                charge = atom.charge
                mass = atom.mass
                if mass == 0.:
                    # need to assign mass
                    # just to trick cpptraj to avoid create empty atom
                    # with mass=0.0
                    mass = 1E-6
            atom = Atom(aname, atype, charge, mass)
            pseudotop.add_atom(atom=atom, resid=resid, resname=resname)

    if is_mdtraj(parm):
        # not sure how to get angles, dihedrals quickly
        if list(parm.bonds):
            pseudotop.add_bonds(np.array([(a.index, b.index) for (a, b) in parm.bonds]))

    elif is_mdanalysis(parm):
        # turn-off. need to check MDAnalysis
        # ack, lots of try and except
        try:
            pseudotop.add_bonds(np.asarray(parm.universe.bonds.to_indices()))
        except TypeError:
            pass

        try:
            pseudotop.add_angles(np.asarray(parm.universe.angles.to_indices()))
        except TypeError:
            pass

        try :
            pseudotop.add_dihedrals(np.asarray(parm.universe.torsions.to_indices()))
        except TypeError:
            pass

        try:
            pseudotop.box = Box(parm.dimensions.astype(np.float64))
        except AttributeError:
            pass
    else:
        # parmed
        if parm.bonds:
            bond_list = [(x.atom1.idx, x.atom2.idx)
                        for x in parm.bonds]
        else:
            bond_list = []

        if parm.angles:
            angle_list = [(x.atom1.idx, x.atom2.idx, x.atom3.idx)
                        for x in parm.angles]
        else:
            angle_list = []

        if parm.dihedrals:
            dihedral_list= [(x.atom1.idx, x.atom2.idx, x.atom3.idx, x.atom4.idx)
                      for x in parm.dihedrals]
        else:
            dihedral_list = []

        if bond_list:
            pseudotop.add_bonds(np.asarray(bond_list))

        if angle_list:
            pseudotop.add_angles(np.asarray(angle_list))

        if dihedral_list:
            pseudotop.add_dihedrals(np.asarray(dihedral_list))

        try:
            pseudotop.box = Box(np.array(parm.box))
        except (ValueError, TypeError):
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
