"""Load mdtraj traj object
"""
from __future__ import absolute_import


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
    from pytraj.api import Trajectory
    from ._load_pseudo_parm import load_pseudo_parm

    if autoconvert:
        unit = 10.
    else:
        unit = 1.

    if not isinstance(m_traj, MDTrajectory):
        raise ValueError("must be mdtraj's Trajectory object")
    else:
        if top is not None:
            pseudotop = top
        else:
            # make Topology, a bit slow
            pseudotop = load_pseudo_parm(m_traj.top)
            if not m_traj.unitcell_lengths is None:
                # convert "nm" to "Angstrom"
                # only check box in 1st frame
                arr = np.append(
                    unit * m_traj.unitcell_lengths[0], m_traj.unitcell_angles[0])
                pseudotop.box = Box(arr.astype(np.float64))

        traj = Trajectory(xyz=m_traj.xyz, top=pseudotop)
        traj.unitcells = np.hstack((unit*m_traj.unitcell_lengths,
                m_traj.unitcell_angles))
        return traj 
