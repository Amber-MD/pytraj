"""perform nucleic acid analysis
"""
from __future__ import absolute_import
import numpy as np
from ._base_result_class import BaseAnalysisResult
from .externals.six import string_types
from ._get_common_objects import _get_topology, _get_data_from_dtype
from ._get_common_objects import _get_reference_from_traj, _get_fiterator
from pytraj.externals.six import iteritems

def _group(self, key):
    # adapted from `toolz` package.
    # see license in $PYTRAJHOME/licenses/externals/toolz.txt
    import collections
    if not callable(key):
        key = getter(key)
    #d = collections.defaultdict(lambda: self.__class__().append)
    d = collections.defaultdict(lambda: [].append)
    for item in self:
        d[key(item)](item)
    rv = {}
    for k, v in iteritems(d):
        rv[k] = v.__self__
    return rv


def nastruct(traj=None,
             ref=0,
             mask="",
             resmap=None,
             hbcut=None,
             frame_indices=None,
             top=None):
    """compute nucleic acid parameters.

    Parameters
    ----------
    traj : Trajectory-like
    ref : {Frame, int}, default 0 (first frame)
    mask : atom mask
    resmap : residue map, example: 'AF2:A'
    hbcut : float
    frame_indices : array-like, default None (all frames)

    Returns
    -------
    out : nupars object. One can assess different values (major groove width, xdips values
    ...) by accessing its attribute. See example below.

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> data = pt.nastruct(traj)
    >>> data.keys()[:5]
    ['buckle', 'minor', 'major', 'xdisp', 'stagger']
    >>> # get minor groove width values for each pairs for each snapshot
    >>> # data.minor is a tuple, first value is a list of basepairs, seconda value is
    >>> # numpy array, shape=(n_frames, n_pairs)

    >>> data.minor
    (['1G16C', '2G15C', '3G14C', '4C13G', '5G12C', '6C11G', '7C10G', '8C9G'], 
     array([[ 13.32927036,  13.403409  ,  13.57159901, ...,  13.26655865,
             13.43054485,  13.4557209 ],
           [ 13.32002068,  13.45918751,  13.63253593, ...,  13.27066231,
             13.42743683,  13.53450871],
           [ 13.34087658,  13.53778553,  13.57062435, ...,  13.29017353,
             13.38542843,  13.46101475]]))

    >>> data.twist
    (['1G16C-2G15C', '2G15C-3G14C', '3G14C-4C13G', '4C13G-5G12C', '5G12C-6C11G', '6C11G-7C10G', '7C10G-8C9G'], 
    array([[ 34.77773666,  33.98158646,  30.18647003, ...,  35.14608765,
             33.9628334 ,  33.13056946],
           [ 33.39176178,  32.68476105,  28.36385536, ...,  36.59774399,
             30.20827484,  26.48732948],
           [ 36.20665359,  32.58955002,  27.47707367, ...,  33.42843246,
             30.90047073,  33.73724365]]))
    """
    from pytraj.datasets.DatasetList import DatasetList as CpptrajDatasetList
    from .actions.CpptrajActions import Action_NAstruct
    from pytraj.array import DataArray

    fi = _get_fiterator(traj, frame_indices)
    _ref = _get_reference_from_traj(traj, ref)
    _top = _get_topology(traj, top)
    _resmap = "resmap " + resmap if resmap is not None else ""
    _hbcut = "hbcut " + str(hbcut) if hbcut is not None else ""

    if not isinstance(mask, string_types):
        # [1, 3, 5] to "@1,3,5
        mask = to_cpptraj_atommask(mask)

    command = " ".join((mask, _resmap, _hbcut))

    act = Action_NAstruct()
    dslist = CpptrajDatasetList()

    act(command, [_ref, fi], dslist=dslist, top=_top)

    # need to update key to avoid duplicate (same key with different
    # aspect)
    dslist_py = []
    for d in dslist:
        # for panda's dataframe
        #d.key = 'nuc_' + d.key + "_" + d.aspect
        # exclude reference value
        dslist_py.append(DataArray(d))
        dslist_py[-1].values = dslist_py[-1].values[1:]
    return nupars(_group(dslist_py, lambda x : x.aspect))


class nupars(object):
    '''class holding data for nucleic acid.
    '''
  
    def __init__(self, adict):
        self._dict = adict

    def __str__(self):
        return '<nupars, keys = %s>' % str(self.keys())

    def __repr__(self):
        return str(self)

    def __getitem__(self, key):
        '''self['minor'], ...
        '''
        return self.__getattr__(key)

    def __getattr__(self, aspect):
        '''self.minor, ...
        '''
        # data is a list of DataArray
        data = self._dict[aspect]
        arr = np.empty((len(data), len(data[0])), dtype='f8')

        keylist = []
        for idx, arr0 in enumerate(data):
            keylist.append(arr0.key)
            arr[idx] = arr0.values
        return keylist, arr.T

    def keys(self):
        return list(self._dict)

    def __dir__(self):
        '''for autocompletion in ipython
        '''
        return self.keys()

    def _summary(self, op, keys=None, indices=None):
        '''
        Parameters
        op : numpy method
        keys: optional
        indices : optional

        Examples
        --------
        self._summary(np.mean, indices=range(2, 8))
        '''
        _keys = keys if keys is not None else self.keys()

        sumdict = {}
        for k in _keys:
            values = self[k][1]
            if indices is None:
                sumdict[k] = op(values, axis=1)
            else:
                sumdict[k] = op(values[indices], axis=1)
        return sumdict
