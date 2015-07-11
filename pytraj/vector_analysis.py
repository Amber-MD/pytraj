from __future__ import print_function, absolute_import

from .externals.six import iteritems

supported_types = [x for x in
                   'mask minimage dipole center corrplane box boxcenter ucellx ucelly ucellz principal'.split()]

template = '''
def vector_%s(traj=None, command="", top=None, *args, **kwd):
    """
    Parameters
    ----------
    traj : Trajectory-like or anything that makes _frame_iter_master(traj) return Frame
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments
    """
    from ._get_common_objects import _get_top, _get_data_from_dtype
    from ._get_common_objects import _get_list_of_commands
    from .datasets.DataSetList import DataSetList
    from .actions.CpptrajActions import Action_Vector
    from .core.ActionList import ActionList

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = ' %s '

    list_of_commands = _get_list_of_commands(command)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add_action(act, _command, _top, dslist=dslist, *args, **kwd)
    actlist.do_actions(traj)
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_types:
    my_func_str = template % (key, key)
    g_dict = globals()
    exec(my_func_str)

del key
