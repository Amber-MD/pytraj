""""""
from __future__ import absolute_import
from ._base_result_class import BaseAnalysisResult
from .externals.six import string_types
from ._get_common_objects import _get_top, _get_data_from_dtype
from ._get_common_objects import _get_reference_from_traj


def nastruct(traj=None,
             mask="",
             resmap=None,
             hbcut=None,
             ref=None,
             top=None,
             dtype='nupars', *args, **kwd):
    """
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.nastruct(traj0)
    >>> pt.nastruct(traj1, resmap='AF2:A')
    """
    from pytraj.datasets.DatasetList import DatasetList as CpptrajDatasetList
    from pytraj.datasetlist import DatasetList as Dataset
    from .actions.CpptrajActions import Action_NAstruct
    from pytraj.array import DataArray

    _ref = _get_reference_from_traj(traj, ref)
    _resmap = "resmap " + resmap if resmap is not None else ""
    _hbcut = "hbcut " + str(hbcut) if hbcut is not None else ""

    if not isinstance(mask, string_types):
        # [1, 3, 5] to "@1,3,5
        mask = to_cpptraj_atommask(mask)

    command = " ".join((mask, _resmap, _hbcut))

    act = Action_NAstruct()
    dslist = CpptrajDatasetList()

    _top = _get_top(traj, top)
    act(command, [_ref, traj], dslist=dslist, top=_top, *args, **kwd)

    # need to update legend to avoid duplicate (same legend with different
    # aspect)
    dslist_py = Dataset()
    for d in dslist:
        # for panda's dataframe
        d.legend = 'nuc_' + d.legend + "_" + d.aspect
        # exclude reference value
        dslist_py.append(DataArray(d))
        dslist_py[-1].values = dslist_py[-1].values[1:]
    if dtype == 'nupars':
        return nupars(_get_data_from_dtype(dslist, dtype='dataset'))
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


class nupars(BaseAnalysisResult):
    '''class holding data for nucleic acid.
    Just use ``mean_and_std`` to get mean and std of each component (major groove, ...)
    '''

    def mean_and_std(self):
        """return a dict
        """
        import numpy as np

        aspects = self._dslist.aspects()

        d = {}
        for idx, ap in enumerate(aspects):
            d0 = self._dslist.grep(ap, mode='aspect')
            d0_mean = d0.mean()[1].mean()
            d0_std = d0.mean()[1].std()
            d[ap] = np.array([d0_mean, d0_std], dtype='f8')
        return d
