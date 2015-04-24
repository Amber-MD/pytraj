from __future__ import absolute_import
from .externals.six import string_types
from ._utils import set_world_silent
from .Topology import Topology
from .DataSetList import DataSetList

def _get_top(traj, top):
    if isinstance(top, string_types):
        _top = Topology(top)
    elif top is None: 
        if hasattr(traj, 'top'):
           _top = traj.top 
        else:
            # list, tuple of traj objects 
            try:
                for tmp in traj:
                    if hasattr(tmp, 'top'):
                        _top = tmp.top 
                        break
            except:
                print("Topology is None")
                _top = None
    else:
        _top = top
    return _top

def calculate(action=None, command=None, traj=None, top=None, 
              dslist=DataSetList(), quick_get=False, **kwd): 
    """ quick way to get data 
    Parameters
    ----------
    action : Action object or str, optional
    command : str, default=None 
        command for specific action. For example, if action=`rmsd`, command might be `@CA`
    traj : Trajectory object (FrameArray, TrajReadOnly, ...) or list, tuple of traj object 
    top : topology 
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
            from pytraj.Trajin_Single import Trajin_Single
            traj = Trajin_Single(traj, _top)
        except:
            raise ValueError("can not load %s" % traj)
    if isinstance(action, string_types): 
        # convert to action 
        act = adict[action] 
    else: 
        act = action 
    act(command, traj, _top, dslist=dslist, quick_get=quick_get, **kwd)

    # make new view for DataSetList
    # for some reasons, if I kept calling `calculate` several times,
    # data will be added to the same dslist
    _dslist = DataSetList()
    _dslist.set_py_free_mem(False)
    for idx, ds in enumerate(dslist):
        if idx >= old_size:
            _dslist.add_existing_set(ds)
    if quick_get:
        return _dslist[-1]
    else:
        return _dslist
