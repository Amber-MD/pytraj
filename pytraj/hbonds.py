from __future__ import absolute_import, print_function

from pytraj.action_dict import ActionDict
from .externals.six import string_types
from pytraj.DataSetList import DataSetList
from pytraj._get_common_objects import _get_data_from_dtype

adict = ActionDict()

def search_hbonds(traj, mask="", dtype='dataset', *args, **kwd):
    """search hbonds for a given mask
    Parameters
    ----------
    traj : {Trajectory-like object, frame_iter object, list of traj}
    mask : str 
        Amber atom mask
    dtype : str {'list', 'pyarray', 'dataset', 'ndarray'}, default='dataset'
    *args, **kwd: optional

    Returns
    -------
    out : DataSetList | pyarray | ndarray | list | dict (depend on 'dtype')

    Examples
    --------
    * The syntax was adapted from http://ambermd.org/doc12/Amber15.pdf (page 575)
    * The explaniation in " " is direct excerpt from this manual

    * "search for all hydrogen bonds within residues 1-22"
        dslist = search_hbonds(traj, ":1-22")
        
    * "search for all hydrogen bonds within residues 1-22, specifying output files"

        dslist = search_hbonds(traj, ":1-22 out nhb.dat avgout avghb.dat", dflist=dflist)
        dflist.write_all_datafile()

    * "search for all hydrogen bonds formed between donors in residue 1 and acceptors in residue 2" 

        dslist = search_hbonds(traj, "donormask :1 acceptormask :2", dtype='ndarray'))
   
    See Also
    --------
    http://ambermd.org/doc12/Amber15.pdf (page 575)
    """
    dslist = DataSetList()
    act = adict['hbond']
    command = "series " + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)

def search_nointramol_hbonds(traj, mask="", dtype='dataset', *args, **kwd):
    """
    See Also
    --------
       search_hbonds
    """
    dslist = DataSetList()
    act = adict['hbond']
    command = "series nointramol" + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)
