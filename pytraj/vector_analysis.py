from __future__ import print_function, absolute_import

from .externals.six import iteritems

supported_types= [x for x in 
'mask minimage dipole center corrplane box boxcenter ucellx ucelly ucellz principal'.split()]

template = '''
def calc_%s(traj=None, command="", top=None, *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like or anything that makes _frame_iter_master(traj) return Frame
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments
    """
    from ._get_common_objects import _get_top, _get_data_from_dtype
    from .DataSetList import DataSetList
    from .actions.Action_Vector import Action_Vector

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = '%s '
    template_command += command 

    act = Action_Vector()
    act(template_command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_types:
    my_func_str = template % (key, key)
    g_dict = globals()
    exec(my_func_str)

del key
