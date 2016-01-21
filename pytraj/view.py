'''module to export pytraj's objects to external objects. Mostly for visulization
'''


def to_chemview(traj):  # pragma: no cover
    '''Generate topology spec for the MolecularViewer from pytraj. (adapted from chemview code)

    Parameters
    ----------
    traj : pytraj.Trajectory or pytraj.TrajectoryIterator
    '''
    import pytraj as pt
    from chemview import TrajectoryViewer

    top = {}
    top['atom_types'] = [a.element[1] for a in traj.topology.atoms]
    top['atom_names'] = [a.name for a in traj.topology.atoms]
    top['bonds'] = traj.topology.bond_indices
    # only calculate dssp for 1st frame?
    top['secondary_structure'] = pt.dssp_allresidues(traj[:1],
                                                     simplified=True)[0]
    top['residue_types'] = [r.name for r in traj.topology.residues]
    top['residue_indices'] = [
        list(range(r.first_atom_index, r.last_atom_index))
        for r in traj.topology.residues
    ]

    return TrajectoryViewer(traj.xyz, top)


def to_nglview(traj, parmfile=None):  # pragma: no cover
    '''convert to nglview object

    Parameters
    ----------
    traj : pytraj.TrajectoryIterator or Trajectory

    Returns
    -------
    nglview.Trajectory

    Notes
    -----
    need ParmEd

    '''
    try:
        from io import StringIO
    except ImportError:
        from cStringIO import StringIO
    import parmed as pmd
    import nglview as nv
    from pytraj.sandbox import to_parmed

    if parmfile is None:
        parm = to_parmed(traj)
    else:
        parm = pmd.load_file(parmfile)
        parm.coordinates = traj[0].xyz
        parm.box = traj[0].box
    x = StringIO()
    parm.write_pdb(x)
    buffer_ = x.getvalue()
    x.close()
    ngl_traj = nv.Trajectory(xyz=traj.xyz, topology=nv.Structure(text=buffer_))
    return nv.TrajectoryViewer(ngl_traj)
