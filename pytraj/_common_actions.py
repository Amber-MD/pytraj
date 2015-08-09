from __future__ import absolute_import
from .externals.six import string_types
from ._set_silent import set_world_silent
from .Topology import Topology
from .datasets.DataSetList import DataSetList as CpptrajDatasetList
from ._get_common_objects import _get_top


def calculate(action=None,
              traj=None,
              command="",
              top=None,
              dslist=CpptrajDatasetList(),
              dtype='dataset',
              quick_get=False, **kwd):
    """ quick way to get data 
    Parameters
    ----------
    action : Action object or str, optional
    command : str, default=None 
        command for specific action. For example, if action=`rmsd`, command might be `@CA`
    traj : Trajectory object (Trajectory, TrajectoryIterator, ...) or list, tuple of traj object 
    top : topology 
    dslist : CpptrajDatasetList
        Hold output
    dtype : str {'pyarray', 'list', 'ndarray', 'dataset'}
    **kwd : additional arguments

    Use `calculate(ahelp=True)` or `calculate(ahelp='action name')` for help 

    Returns
    -------
    DatSet object

    >>> from pytraj import calculate
    >>> from pytraj import CpptrajDatasetList 
    >>> dslist = CpptrajDatasetList()
    >>> d0 = calculate('distance', ":2@CA :4@CA", traj, dslist=dslist)
    >>> # d0 == dslist[-1]

    """
    from pytraj import ActionDict
    adict = ActionDict()

    not_return_list = ['action_autoimage', 'action_translate', 'action_rotate',
                       'action_scale', 'action_strip']

    need_print_output = False
    if kwd:
        if 'print_output' in kwd.keys() and kwd['print_output'] == True:
            need_print_output = True
            # need to remove "print_output" since Action does not have this
            # keyword
        else:
            need_print_output = False
        if 'print_output' in kwd.keys():
            kwd.pop('print_output', None)

    if dslist is None:
        dslist = CpptrajDatasetList()
    elif not isinstance(dslist, CpptrajDatasetList):
        raise NotImplementedError(
            "must have None or CpptrajDatasetList object")

    old_size = dslist.size

    if traj is None or isinstance(traj, string_types):
        raise ValueError("must have trajectory object")

    try:
        _top = _get_top(traj, top)
    except:
        raise ValueError("can not get Topology")

    if isinstance(action, string_types):
        # convert to action
        act = adict[action]
    else:
        act = action

    act_name = act.__class__.__name__.lower()

    act(command, traj, _top, dslist=dslist, quick_get=quick_get, **kwd)
    if need_print_output:
        act.print_output()

    # make new view for CpptrajDatasetList
    # for some reasons, if I kept calling `calculate` several times,
    # data will be added to the same dslist
    _dslist = CpptrajDatasetList()
    _dslist.set_py_free_mem(False)
    for idx, ds in enumerate(dslist):
        if idx >= old_size:
            _dslist.add_existing_set(ds)
    if quick_get:
        _d = dslist[-1]
        if dtype == 'pyarray':
            return _d.to_pyarray()
        elif dtype == 'list':
            return _d.tolist()
        elif dtype == 'ndarray':
            return _d.to_ndarray()
        elif dtype == 'dataset':
            return _d
        else:
            raise NotImplementedError(dtype)
    else:
        if act_name not in not_return_list:
            return _dslist
        else:
            return None
