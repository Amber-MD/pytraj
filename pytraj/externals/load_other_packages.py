from __future__ import absolute_import
from ..utils.get_common_objects import _load_Topology
from ..utils.context import tempfolder
from .six import string_types


def load_parmed(parm, traj=True, **kwd):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ----------
    parm : ParmEd's Structure object
    traj: bool, default True
        if True, return pytraj.Trajectory
        if False, return Topology

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> traj = pt.load_parmed(p)
    """
    from parmed.amber import AmberParm

    with tempfolder():
        if isinstance(parm, AmberParm):
            fname = 'tmp.parm7'
        else:
            fname = 'tmp.psf'
        parm.save(fname)
        top = _load_Topology(fname)
    if traj:
        from pytraj import Trajectory
        coords = parm.get_coordinates()
        return Trajectory(xyz=coords, top=top)
    else:
        return top


def load_mdtraj(m_traj, autoconvert=False, top=None):
    """load_mdtraj traj object

    Parameters
    ----------
    m_traj : Trajectory object from mdtraj
    autoconvert : bool, default=False
        convert from "nm" (mdtraj )to "Angstrom" (pytraj)
    top : pytraj.Topology, default None
        if `top` is provided, the converting will be much faster
    """
    import numpy as np
    from mdtraj import Trajectory as MDTrajectory
    from pytraj.trajectory import Trajectory

    if autoconvert:
        unit = 10.
    else:
        unit = 1.

    if not isinstance(m_traj, MDTrajectory):
        raise ValueError("must be mdtraj's Trajectory object")
    else:
        if isinstance(top, string_types):
            pseudotop = _load_Topology(top)
        else:
            pseudotop = top
        if pseudotop is None:
            raise ValueError("need Topology or pdb, mol2, ... files")
        traj = Trajectory(xyz=m_traj.xyz.astype('f8'), top=pseudotop)

        if m_traj.unitcell_lengths is not None and m_traj.unitcell_angles is not None:
            traj.unitcells = np.hstack((unit * m_traj.unitcell_lengths,
                                        m_traj.unitcell_angles)).astype('f8')
        return traj


def load_MDAnalysis(universe, top=None):
    """load MDAnalysis' Universe object to pytra's traj object

    Notes
    -----
    All coords will be loaded

    See Also
    --------
    load_MDAnalysisIterator
    """

    from MDAnalysis import Universe
    from pytraj.Trajectory import Trajectory
    from ..Frame import Frame

    # don't import here since we import load_pseudo_parm in
    # TrajectoryMDAnalysisIterator

    # MDAnalysis needs numpy. So we always have numpy when using this
    if not isinstance(universe, Universe):
        raise ValueError("must be a Universe")

    # creat pseudotop
    if top is None:
        raise ValueError("need a Topology or pdb/mol2/...")
    else:
        pseudotop = top

    # creat atom group
    ag = universe.atoms

    farray = Trajectory()
    farray.top = pseudotop
    for _ in universe.trajectory:
        frame = Frame(farray.top.n_atoms)
        # set box for each Frame
        frame.boxview[:] = farray.top.box[:]
        # load xyz coords, let numpy do automatically casting
        frame.xyz[:] = ag.positions
        # we don't need to make copy=True since we already created
        # frame and `farray` can 'keep' it
        farray.append(frame, copy=False)
    return farray
