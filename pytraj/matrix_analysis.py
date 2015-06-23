from __future__ import print_function, absolute_import
from .externals.six import iteritems


__all__ = ['distance_matrix', 'correlation_matrix', 'coord_covariance_matrix',
           'mw_covariance_matrix', 'distcovar_matrix', 'idea_matrix',
           ]

default_key_dict = {'distance_matrix' : 'dist',
        'correlation_matrix' : 'correl',
        'coord_covariance_matrix' : 'covar',
        'mw_covariance_matrix' : 'mwcovar',
        'distcovar_matrix' : 'distcovar',
        'idea_matrix' : 'idea'}

__cpptrajdoc__ = """
    cpptraj manual
    --------------
          [out <filename>] [start <start>] [stop <stop>] [offset <offset>]
          [name <name>] [ byatom | byres [mass] | bymask [mass] ]
          [ ired [order <#>] ]
          [ {distcovar | idea} <mask1> ]
          [ {dist | correl | covar | mwcovar} <mask1> [<mask2>] ]
          [ dihcovar dihedrals <dataset arg> ]
    Calculate a matrix of the specified type from input coordinates.
      dist: Distance matrix (default).
      correl: Correlation matrix (aka dynamic cross correlation).
      covar: Coordinate covariance matrix.
      mwcovar: Mass-weighted coordinate covariance matrix.
      distcovar: Distance covariance matrix.
      idea: Isotropically Distributed Ensemble Analysis matrix.
      ired: Isotropic Reorientational Eigenmode Dynamics matrix.
      dihcovar: Dihedral covariance matrix.
"""

template = '''
def %s(traj=None, command="", top=None, dtype='ndarray', *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like or anything that makes _frame_iter_master(traj) return Frame
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments

    Examples
    --------
    >>> from pytraj import matrix_analysis as ma
    >>> from pytraj import io
    >>> traj = io.load_sample_data('tz2')
    >>> ma.distcovar_matrix(traj, '@CA')
    >>> dslist = ma.mw_covariance_matrix(traj, '@CA')
    >>> #print (dslist) 

    See Also
    --------
    default_key_dict
    _shared_methods._frame_iter_master

    cpptraj compat mode
    -------------------
    {'distance_matrix' : 'dist',
    'correlation_matrix' : 'correl',
    'coord_covariance_matrix' : 'covar',
    'mw_covariance_matrix' : 'mwcovar',
    'distcovar_matrix' : 'distcovar',
    'idea_matrix' : 'idea'}
    """
    from .actions.CpptrajActions import Action_Matrix
    from ._get_common_objects import _get_top, _get_data_from_dtype
    from .datasets.DataSetList import DataSetList

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = '%s '
    template_command += command 

    act = Action_Matrix()
    act(template_command, traj, top=_top, dslist=dslist, *args, **kwd)
    # need to call `print_output` so cpptraj can normalize some data
    # check cpptraj's code
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for k, v in iteritems(default_key_dict):
    my_func_str = template % (k, v)
    g_dict = globals()
    exec(my_func_str)
    g_dict[k].__doc__ += __cpptrajdoc__

del k, v
