from __future__ import absolute_import, print_function

from . import adict
from .externals.six import string_types

def search_hbonds(traj, mask="", *args, **kwd):
    """return a DataSetList object storing all pair of h-bond
    Parameters
    ----------
    mask : str (Amber atom mask)
    traj : {Trajectory-like object, frame_iter object, list of traj}
    top : optional, Topology file

    *args, **kwd: optional

    Returns
    ------
    DataSetList object
    """
    act = adict['hbond']
    command = "series " + mask
    return act(command, traj, *args, **kwd)
