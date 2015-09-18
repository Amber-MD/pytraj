from __future__ import print_function, absolute_import

from .externals.six import iteritems

supported_types = [
    x for x in
    'mask minimage dipole center corrplane box boxcenter ucellx ucelly ucellz principal'.split(
    )
]

template = '''
def vector_%s(traj=None, command="", frame_indices=None, dtype='ndarray', top=None):
    """
    Parameters
    ----------
    traj : Trajectory-like
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments
    """
    from ._get_common_objects import _get_top, _get_data_from_dtype, _get_fiterator
    from ._get_common_objects import _get_list_of_commands
    from .datasets.DatasetList import DatasetList as CpptrajDatasetList
    from .actions.CpptrajActions import Action_Vector
    from .core.ActionList import ActionList

    fi = _get_fiterator(traj, frame_indices)
    _top = _get_top(fi, top)
    dslist = CpptrajDatasetList()
    template_command = ' %s '

    list_of_commands = _get_list_of_commands(command)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add_action(act, _command, _top, dslist=dslist)
    actlist.do_actions(fi)
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_types:
    my_func_str = template % (key, key)
    g_dict = globals()
    exec(my_func_str)

del key
