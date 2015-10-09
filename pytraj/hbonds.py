from __future__ import absolute_import, print_function, division

from .action_dict import ActionDict
from .externals.six import string_types
from .datasets.DatasetList import DatasetList
from ._get_common_objects import _get_data_from_dtype, _get_topology
from ._base_result_class import BaseAnalysisResult

adict = ActionDict()

__all__ = ['DatasetHBond', 'search_hbonds', 'search_nointramol_hbonds',
           'search_hbonds_noseries']


def _to_amber_mask(txt):
    import re
    """Convert something like 'ASP_16@OD1-ARG_18@N-H to ':16@OD1 :18@H'
    """

    if isinstance(txt, string_types):
        _txt = [txt, ]
    elif isinstance(txt, (list, tuple)):
        _txt = txt[:]
    else:
        raise NotImplementedError()

    for mask in _txt:
        mask = mask.replace("_", ":")
        yield " ".join(re.findall(r"(:\d+@\w+)", mask))


class DatasetHBond(BaseAnalysisResult):
    """Hold data for hbond analysis

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_pdb_rcsb("1l2y")
    >>> h = pt.hbonds.DatasetHBond(traj.search_hbonds())
    >>> h
    <pytraj.hbonds.DatasetHBond
    donor_aceptor pairs : 31>
    """

    def __str__(self):
        root_msg = "<pytraj.hbonds.DatasetHBond"
        more_info = "donor_aceptor pairs : %s>" % len(self.donor_aceptor)
        return root_msg + "\n" + more_info

    def __repr__(self):
        return self.__str__()

    @property
    def donor_aceptor(self):
        return self.dataset.grep(["solventhb", "solutehb"],
                                 mode='aspect').keys()

    def to_amber_mask(self):
        """convert donor_aceptor pair mask to amber mask to calculate
        distance (for example: 'ARG16_O-TRP6_NE1-HE1' will be ':16@O :6@HE1')
        """
        raise NotImplementedError("not yet")


def _update_key_hbond(_dslist):

    # SER_20@O-SER_20@OG-HG --> SER20_O-SER20_OG-HG
    for d0 in _dslist:
        d0.key = d0.key.replace("_", "")
        d0.key = d0.key.replace("@", "_")

    for d0 in _dslist:
        if d0.key == 'HB00000[UU]':
            d0.key = 'total_solute_hbonds'


def search_hbonds_noseries(traj,
                           mask="",
                           dtype='dataset',
                           update_key=True, *args, **kwd):
    """search hbonds for a given mask

    Parameters
    ----------
    traj : {Trajectory-like object, frame_iter object, list of traj}
    mask : str 
        Amber atom mask
    dtype : str {'dataset', 'ndarray'}, default='dataset'
    *args, **kwd: optional

    Returns
    -------
    out : pytraj.DatasetList (similiar to Python list with labeled array)

    See Also
    --------
    http://ambermd.org/doc12/Amber15.pdf (page 575)
    """

    dslist = DatasetList()
    act = adict['hbond']

    command = mask
    if "series" in command:
        raise ValueError("don't accept key `series`")
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()

    if update_key:
        _update_key_hbond(dslist)

    if dtype == 'dataframe':
        # return DataFrame.T to have better visual effect
        return dslist.to_dataframe().T
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def search_hbonds(traj,
                  mask="",
                  solvent_donor=None,
                  solvent_acceptor=None,
                  distance=3.0,
                  angle=135.,
                  dtype='hbond',
                  image=False,
                  more_options='',
                  top=None):
    """search hbonds for a given mask. 

    Parameters
    ----------
    traj : Trajectory-like
    mask : str 
        Amber atom mask

    Returns
    -------
    out : pytraj.DatasetList, which is similiar to Python list with labeled data.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> pt.search_hbonds(traj, ':5,8')
    <pytraj.DatasetList with 3 datasets>
    total_solute_hbonds
    [2 2 0 1 1 0 2 1 1 1]
    LYS8_O-GLU5_N-H
    [1 1 0 1 1 0 1 1 1 1]
    GLU5_O-LYS8_N-H
    [1 1 0 0 0 0 1 0 0 0]
    """
    from pytraj.actions.CpptrajActions import Action_Hbond
    dslist = DatasetList()
    act = Action_Hbond()

    _top = _get_topology(traj, top)

    s_donor = "solventdonor " + str(solvent_donor) if solvent_donor else ""
    s_acceptor = "solventacceptor " + \
        str(solvent_acceptor) if solvent_acceptor else ""
    _dist = 'dist ' + str(distance)
    _angle = 'angle ' + str(angle)
    _image = 'image' if image else ''
    _options = more_options

    command = " ".join(
        ("series", mask, s_donor, s_acceptor, _dist, _angle, _image, _options))
    act(command, traj, top=_top, dslist=dslist)
    act.print_output()

    old_keys = dslist.keys()
    _update_key_hbond(dslist)
    if dtype == 'dataframe':
        # return DataFrame.T to have better visual effect
        return dslist.to_dataframe().T
    elif dtype == 'hbond':
        dslist_new = _get_data_from_dtype(dslist, dtype='dataset')
        hdata = DatasetHBond(dslist_new)
        hdata._old_keys = old_keys
        return hdata
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def search_nointramol_hbonds(traj,
                             mask="solventacceptor :WAT@O solventdonor :WAT",
                             dtype='dataset',
                             update_key=True, *args, **kwd):
    """
    Search hbonds between solute and solvent, ignoring intra-hbond

    Parameters
    ----------
    traj : Trajectory-like or any iterable object that _frame_iter_mater return a Frame
    mask : str, default "solventacceptor :WAT@O solventdonor :WAT"
        cpptraj command
    dtype : str, default 'dataset'
    *args, **kwd: optional

    Examples
    --------
    >>> import pytraj as pt
    >>> pt.search_nointramol_hbonds(traj)
    >>> pt.search_nointramol_hbonds([traj, traj2], top=traj.top)

    See Also
    --------
       search_hbonds
    """
    from pytraj.actions.CpptrajActions import Action_Hbond
    dslist = DatasetList()
    act = Action_Hbond()
    command = "series nointramol " + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()

    if update_key:
        _update_key_hbond(dslist)
    if dtype == 'hbond_class':
        return HbondAnalysisResult(dslist)
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)
