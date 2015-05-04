from __future__ import absolute_import
from .externals.six import string_types
from ._set_silent import set_world_silent
from .Topology import Topology
from .DataSetList import DataSetList
from ._get_common_objects import _get_top

def calculate(action=None, traj=None, command="", top=None, 
              dslist=DataSetList(), dtype='dataset', quick_get=False, **kwd): 
    """ quick way to get data 
    Parameters
    ----------
    action : Action object or str, optional
    command : str, default=None 
        command for specific action. For example, if action=`rmsd`, command might be `@CA`
    traj : Trajectory object (FrameArray, TrajReadOnly, ...) or list, tuple of traj object 
    top : topology 
    dslist : DataSetList
        Hold output
    dtype : str {'pyarray', 'list', 'ndarray', 'dataset'}
    **kwd : additional arguments
 
    Use `calculate(ahelp=True)` or `calculate(ahelp='action name')` for help 

    Returns
    -------
    DatSet object

    >>> from pytraj import calculate
    >>> from pytraj import DataSetList 
    >>> dslist = DataSetList()
    >>> d0 = calculate('distance', ":2@CA :4@CA", traj, dslist=dslist)
    >>> # d0 == dslist[-1]
 
    """ 
    from pytraj import ActionDict
    adict = ActionDict()

    need_print_output = False
    if kwd:
        if 'print_output' in kwd.keys() and kwd['print_output'] == True:
            need_print_output = True
            # need to remove "print_output" since Action does not have this keyword
        else:
            need_print_output = False
        if 'print_output' in kwd.keys():
            kwd.pop('print_output', None)

    if dslist is None:
        dslist = DataSetList()
    elif not isinstance(dslist, DataSetList):
        raise NotImplementedError("must have None or DataSetList object")

    old_size = dslist.size

    _top = _get_top(traj, top)

    if traj is None: 
        raise ValueError("must have trajectory object") 
    elif isinstance(traj, string_types):
        try:
            from pytraj.trajs.Trajin_Single import Trajin_Single
            traj = Trajin_Single(traj, _top)
        except:
            raise ValueError("can not load %s" % traj)
    if isinstance(action, string_types): 
        # convert to action 
        act = adict[action] 
    else: 
        act = action 
    act(command, traj, _top, dslist=dslist, quick_get=quick_get, **kwd)
    if need_print_output:
        act.print_output()

    # make new view for DataSetList
    # for some reasons, if I kept calling `calculate` several times,
    # data will be added to the same dslist
    _dslist = DataSetList()
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
        return _dslist
