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
    from pytraj.Topology import Topology
    from pytraj.compat import string_types

    if autoconvert:
        unit = 10.
    else:
        unit = 1.

    if not isinstance(m_traj, MDTrajectory):
        raise ValueError("must be mdtraj's Trajectory object")
    else:
        if isinstance(top, string_types):
            pseudotop = Topology(top)
        else:
            pseudotop = top
        if pseudotop is None:
            raise ValueError("need Topology or pdb/mol2/... files")
        traj = Trajectory(xyz=m_traj.xyz, top=pseudotop)
        traj.unitcells = np.hstack((unit*m_traj.unitcell_lengths,
                m_traj.unitcell_angles))
        return traj 
