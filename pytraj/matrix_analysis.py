from __future__ import print_function, absolute_import
from .externals.six import iteritems

__all__ = ['distance_matrix',
           'correlation_matrix',
           'coord_covariance_matrix',
           'mw_covariance_matrix',
           'distcovar_matrix',
           'idea_matrix', ]

default_key_dict = {
    'distance_matrix': 'dist',
    'correlation_matrix': 'correl',
    'coord_covariance_matrix': 'covar',
    'mw_covariance_matrix': 'mwcovar',
    'distcovar_matrix': 'distcovar',
    'idea_matrix': 'idea'
}

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
def %s(traj=None, command="", top=None, dtype='ndarray', mat_type='full', *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like
    command : cpptraj command
    top : Topology, optional, default None
    mat_type : str, {'full', 'half', 'cpptraj'}, default 'full'
        if 'full': 2D full matrix
        if 'half': triangular matrix
        if 'cpptraj': 1D array
    *args, **kwd: more arguments

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
    if dtype == 'ndarray':
        if mat_type == 'full':
            return dslist[0].values
        elif mat_type == 'half':
            return dslist[0].to_half_matrix()
        elif mat_type == 'cpptraj':
            return dslist[0].to_cpptraj_sparse_matrix()
        else:
            raise ValueError()
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)
'''

for k, v in iteritems(default_key_dict):
    my_func_str = template % (k, v)
    g_dict = globals()
    exec(my_func_str)
    g_dict[k].__doc__ += __cpptrajdoc__

del k, v
