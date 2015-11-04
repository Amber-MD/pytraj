from __future__ import absolute_import, print_function, division

from .actions import CpptrajActions
from .action_dict import ActionDict
from .externals.six import string_types
from .datasets import CpptrajDatasetList
from .decorators import _register_pmap
from ._get_common_objects import _get_data_from_dtype, _get_topology
from .base_holder import BaseDataHolder
from ._shared_methods import iterframe_master

__all__ = ['DatasetHBond', 'search_hbonds']


def to_amber_mask(txtlist):
    """Convert hbond lables to amber mask, example 'ASP_16@OD1-ARG_18@N-H to ':16@OD1 :18@H'.
    This converter is good to measure the hbond distance over time.

    >>> list(to_amber_mask(['ASP_16@OD1-ARG_18@N-H',]))
    [':16@OD1 :18@H']
    """
    # TODO: it's better to have cpptraj printing the atom indices
    import re
    # make a copy
    _txt = txtlist[:]

    for mask in _txt:
        mask = mask.replace("_", " ").replace("-", " ").split()
        second_res = mask[3].split('@')[0]
        yield ''.join((':', mask[1], ' :', "".join((second_res, '@', mask[4]))))


class DatasetHBond(BaseDataHolder):
    """Hold data for hbond analysis
    """

    def __str__(self):
        root_msg = "<pytraj.hbonds.DatasetHBond"
        more_info = "donor_aceptor pairs : %s>" % len(self.donor_aceptor)
        return root_msg + "\n" + more_info

    def __repr__(self):
        return str(self)

    @property
    def donor_aceptor(self):
        return self.data.grep(["solventhb", "solutehb"], mode='aspect').keys()

    def _amber_mask(self):
        return list(to_amber_mask(self._old_keys[1:]))


def _update_key_hbond(_dslist):

    # SER_20@O-SER_20@OG-HG --> SER20_O-SER20_OG-HG
    for d0 in _dslist:
        d0.key = d0.key.replace("_", "")
        d0.key = d0.key.replace("@", "_")

    for d0 in _dslist:
        if d0.key == 'HB00000[UU]':
            d0.key = 'total_solute_hbonds'


@_register_pmap
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
    """Searching for Hbond donors/acceptors in region specified by ``mask``

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
    >>> data = pt.search_hbonds(traj, ':5,8')
    >>> data
    <pytraj.hbonds.DatasetHBond
    donor_aceptor pairs : 2>
    >>> data.donor_aceptor
    ['LYS8_O-GLU5_N-H', 'GLU5_O-LYS8_N-H']
    >>> data.values
    array([[2, 2, 0, ..., 1, 1, 1],
           [1, 1, 0, ..., 1, 1, 1],
           [1, 1, 0, ..., 0, 0, 0]], dtype=int32)
    """
    dslist = CpptrajDatasetList()
    act = CpptrajActions.Action_Hbond()

    _top = _get_topology(traj, top)

    s_donor = "solventdonor " + str(solvent_donor) if solvent_donor else ""
    s_acceptor = "solventacceptor " + \
        str(solvent_acceptor) if solvent_acceptor else ""
    _dist = 'dist ' + str(distance)
    _angle = 'angle ' + str(angle)
    _image = 'image' if image else ''
   # TODO: use series=True in args?
    _series = 'series'
    _options = more_options

    command = " ".join(
        (_series, mask, s_donor, s_acceptor, _dist, _angle, _image, _options))

    # need to get correct frame number
    act.read_input(command, top=_top, dslist=dslist)
    act.process(_top)

    for idx, frame in enumerate(iterframe_master(traj)):
        act.do_action(frame, idx=idx)

    act.post_process()

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
