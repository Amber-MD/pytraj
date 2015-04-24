"""Load mdtraj traj object
"""
from __future__ import absolute_import
from pytraj.utils import has_, require
from pytraj.FrameArray import FrameArray
from ._load_pseudo_parm import load_pseudo_parm

def load_mdtraj(m_traj):
    """load_mdtraj traj object

    Parameters
    ---------
    m_traj : Trajectory object from mdtraj 
    """
    if not has_("mdtraj"):
        # we dont need checking `numpy` since mdtraj needs numpy 
        require("mdtraj")
    else:
        from mdtraj import Trajectory
        if not isinstance(m_traj, Trajectory):
            raise PyTrajRequireObject("Trajectory")
        else:
            pseudotop = load_pseudo_parm(m_traj.top)
            return FrameArray(m_traj.xyz, pseudotop)
