"""Load mdtraj traj object
"""
from __future__ import absolute_import
from pytraj.utils import has_, require, _import_numpy
from pytraj.Trajectory import Trajectory
from ._load_pseudo_parm import load_pseudo_parm
from ..Frame import Frame

_, np = _import_numpy()

def load_mdtraj(m_traj, autoconvert=True, top=None):
    """load_mdtraj traj object

    Parameters
    ---------
    m_traj : Trajectory object from mdtraj 
    autoconvert : bool, default=True
        convert from "nm" (mdtraj )to "Angstrom" (pytraj)
    """
    from pytraj.core import Box
    if autoconvert:
        unit = 10.
    else:
        unit = 1.
    if not has_("mdtraj"):
        # we dont need checking `numpy` since mdtraj needs numpy 
        require("mdtraj")
    else:
        from mdtraj import Trajectory as MDTrajectory
        if not isinstance(m_traj, MDTrajectory):
            raise PyTrajRequireObject("Trajectory")
        else:
            if top is not None:
                print ("test")
                pseudotop = top
            else:
                pseudotop = load_pseudo_parm(m_traj.top)
                if not m_traj.unitcell_lengths is None:
                    # convert "nm" to "Angstrom"
                    # only check box in 1st frame
                    arr = np.append(unit*m_traj.unitcell_lengths[0], m_traj.unitcell_angles[0])
                    pseudotop.box = Box(arr.astype(np.float64))

            farray = Trajectory()
            farray.top = pseudotop
            for arr0 in m_traj.xyz:
                frame = Frame(m_traj.n_atoms)
                # convert "nm" to "Angstrom"
                # update xyz for frame
                # make sure to use `float64`
                # TODO: more type-checking
                frame[:] = unit * arr0.astype(np.float64)
                # set box for each Frame
                frame.box = farray.top.box.copy()
                farray.append(frame)
            return farray
