from __future__ import print_function, absolute_import

from functools import partial
from .actions.Action_Matrix import Action_Matrix
from ._get_common_objects import _get_top, _get_data_from_dtype
from .externals.six import iteritems
from .DataSetList import DataSetList

# TODO: some assertions failed. FIXME

__all__ = ['distance_matrix', 'correlation_matrix', 'coord_covariance_matrix',
           'mw_covariance_matrix', 'distcovar_matrix', 'idea_matrix',
           ]

default_key_dict = {'distance_matrix' : 'dist',
        'correlation_matrix' : 'correl',
        'coord_covariance_matrix' : 'covar',
        'mw_covariance_matrix' : 'mwcovar',
        'distcovar_matrix' : 'distcovar',
        'idea_matrix' : 'idea'}

template = """
def %s (traj=None, command="", top=None, *args, **kwd):

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = '%s '
    template_command += command 

    act = Action_Matrix()
    act(template_command, traj, top=_top, dslist=dslist, *args, **kwd)

    return _get_data_from_dtype(dslist, dtype=dtype)
"""

for k, v in iteritems(default_key_dict):
    my_func_str = template % (k, v)
    exec(my_func_str)

del k, v
