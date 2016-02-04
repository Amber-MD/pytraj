from __future__ import absolute_import
from ..get_common_objects import _load_Topology
from ..utils.context import tempfolder
from .six import string_types


def load_ParmEd(parmed_obj, as_traj=False, **kwd):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ---------
    parmed_obj : ParmEd's Structure object
    as_traj: bool, default False
        if True, return pytraj.trajectory.Trajectory
        if False, return Topology

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> top = pt.load_ParmEd(p)
    """
    import parmed as pmd
    from parmed.amber import AmberParm

    if isinstance(parmed_obj, string_types):
        # reserve **kwd for `structure=True`
        parmed_obj = pmd.load_file(parmed_obj, **kwd)
    # faster
    with tempfolder():
        if isinstance(parmed_obj, AmberParm):
            fname = 'tmp.parm7'
        else:
            fname = 'tmppdb.psf'
        parmed_obj.save(fname)
        top = _load_Topology(fname)
    if as_traj:
        from pytraj import Trajectory
        coords = parmed_obj.coordinates
        coords = coords.reshape(1, *coords.shape)
        return Trajectory(xyz=coords, top=top)
    else:
        return top


def _load_parmed(parm_name):
    from parmed import load_file
    return load_file(parm_name)


def to_ParmEd(pytraj_top):
    # TODO: exten to gromacs, charmm too
    # need to change extension
    """convert to ParmEd object"""
    from pytraj.utils.context import tempfolder
    from pytraj.Topology import ParmFile
    import parmed as chem

    # I am not a fan of saving/loading again but this might be best choice
    with tempfolder():
        fname = "tmp_pytrajtop.prmtop"
        ParmFile().writeparm(pytraj_top, fname, format="")
        return chem.load_file(fname)


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
    from pytraj.compat import string_types

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


def load_MDAnalysis(its_obj, top=None):
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
    if not isinstance(its_obj, Universe):
        raise ValueError("must be a Universe")

    # creat pseudotop
    if top is None:
        raise ValueError("need a Topology or pdb/mol2/...")
    else:
        pseudotop = top

    # creat atom group
    ag = its_obj.atoms

    farray = Trajectory()
    farray.top = pseudotop
    for _ in its_obj.trajectory:
        frame = Frame(farray.top.n_atoms)
        # set box for each Frame
        frame.boxview[:] = farray.top.box[:]
        # load xyz coords, let numpy do automatically casting
        frame.xyz[:] = ag.positions
        # we don't need to make copy=True since we already created
        # frame and `farray` can 'keep' it
        farray.append(frame, copy=False)
    return farray
