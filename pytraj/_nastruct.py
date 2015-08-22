""""""
from __future__ import absolute_import
from ._base_result_class import BaseAnalysisResult
from .externals.six import string_types
from ._get_common_objects import _get_top, _get_data_from_dtype
from ._get_common_objects import _get_reference_from_traj


class nupars(BaseAnalysisResult):
    def summary(self):
        """return numpy 2D-array, shape=(n_aspects, 3).
        Each sub-array is [aspect_name, mean, std]

        Examples
        --------
        >>> d.summary()
        [['shear' '-0.034849288128316405' '0.333683347834']
         ['minor' '13.385222816467286' '0.0661943612623']
         ['slide' '-1.6436443726221721' '0.226245717793']
         ..., 
         ['twist' '30.922476874457466' '1.85829981963']
         ['htwist' '32.23108376397027' '1.64625014029']
         ['hb' '2.5300000000000002' '0.483838816136']]
        """
        import numpy as np

        aspects = self._dslist.aspects()

        out = np.empty((len(aspects), 3), dtype='object')

        for idx, ap in enumerate(aspects):
            d0 = self._dslist.grep(ap, mode='aspect')
            d0_mean = d0.mean()[1].mean()
            d0_std = d0.mean()[1].std()
            out[idx] = np.array([ap, d0_mean, d0_std])
        return out


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
    from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
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
