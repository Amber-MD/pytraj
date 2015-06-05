from __future__ import print_function, absolute_import

from .actions.Action_MultiDihedral import Action_MultiDihedral
from ._get_common_objects import _get_top, _get_data_from_dtype
from .externals.six import iteritems
from .DataSetList import DataSetList


__all__ = []

supported_dihedral_types = [x for x in 
'multidihedral phi psi chip omega alpha beta gamma delta epsilon zeta nu1 nu2 chin'.split()]

template = '''
def calc_%s(traj=None, command="", top=None, *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like or anything that makes _frame_iter_master(traj) return Frame
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments

    Examples
    --------
    >>> import pytraj.dihedral_analysis as da
    >>> da.calc_phi(traj)
    >>> da.calc_psi(traj, "resrange 3-10")
    >>> da.calc_chip(traj, "resrange 3-10")
    >>> da.calc_chip(traj, "resrange 3-10", dtype='dict')
    >>> da.calc_multidihedral(traj, "resrange 3-10")
    >>> # assert
    >>> from pytraj import common_actions as pyca
    >>> phi0 = pyca.calc_multidihedral(traj, "phi", dtype='dataset')
    >>> phi1 = da.calc_phi(traj, dtype='dataset')
    >>> from pytraj.testing import aa_eq
    >>> for key in phi0.keys():
    >>>     aa_eq(phi0[key], phi1[key])
    >>> print ("OK")

    See Also
    --------
    """

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = '%s '
    template_command += command 

    act = Action_MultiDihedral()
    act(template_command, traj, top=_top, dslist=dslist, *args, **kwd)
    # need to call `print_output` so cpptraj can normalize some data
    # check cpptraj's code
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_dihedral_types:
    if key != 'multidihedral':
        my_func_str = template % (key, key)
    else:
        my_func_str = template % (key, " ")
    g_dict = globals()
    exec(my_func_str)

del key
