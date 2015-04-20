from __future__ import absolute_import, print_function

from pytraj.action_dict import ActionDict
from .externals.six import string_types
from pytraj.DataSetList import DataSetList

adict = ActionDict()

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
    dslist = DataSetList()
    act = adict['hbond']
    command = "series " + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()
    return dslist

def search_nointramol_hbonds(traj, mask="", *args, **kwd):
    dslist = DataSetList()
    act = adict['hbond']
    command = "series nointramol" + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()
    return dslist
