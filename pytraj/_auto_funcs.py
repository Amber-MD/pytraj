from __future__ import print_function, absolute_import

from ._get_common_objects import _get_top, _get_data_from_dtype
from .externals.six import iteritems
from .DataSetList import DataSetList
from pytraj import adict

# auto-create `calc_` methods


supported_keys = adict.keys()

template = '''
def calc_%s(traj=None, command="", top=None, *args, **kwd):
    from pytraj.action_dict import ActionDict
    from pytraj import DataSetList
    from ._get_common_objects import _get_top, _get_data_from_dtype

    if 'dtype' in kwd.keys():
        dtype = kwd['dtype']
        del kwd['dtype']
    else:
        dtype = None

    _top = _get_top(traj, top)
    dslist = DataSetList()
    template_command = " "
    template_command += command 

    act = ActionDict()["%s"]
    act(template_command, traj, top=_top, dslist=dslist, *args, **kwd)
    # need to call `print_output` so cpptraj can normalize some data
    # check cpptraj's code
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_keys:
    act_name = key
    my_func_str = template % (key, act_name)
    g_dict = globals()
    exec(my_func_str)

del key
