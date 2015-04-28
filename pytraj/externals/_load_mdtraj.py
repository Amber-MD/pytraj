"""Load mdtraj traj object
"""
from __future__ import absolute_import
from pytraj.utils import has_, require, _import_numpy
from pytraj.FrameArray import FrameArray
from ._load_pseudo_parm import load_pseudo_parm
_, np = _import_numpy()

def load_mdtraj(m_traj):
    """load_mdtraj traj object

    Parameters
    ---------
    m_traj : Trajectory object from mdtraj 
    """
    from pytraj.core import Box
    if not has_("mdtraj"):
        # we dont need checking `numpy` since mdtraj needs numpy 
        require("mdtraj")
    else:
        from mdtraj import Trajectory
        if not isinstance(m_traj, Trajectory):
            raise PyTrajRequireObject("Trajectory")
        else:
            pseudotop = load_pseudo_parm(m_traj.top)
            # convert "nm" to "Angstrom"
            fa = FrameArray(10*m_traj.xyz, pseudotop)
            arr = 10 *np.append(m_traj.unitcell_lengths, m_traj.unitcell_angles)
            fa.top.box = Box(arr.astype(np.float64))
            return fa
