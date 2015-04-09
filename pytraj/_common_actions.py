from __future__ import absolute_import
from .externals.six import string_types
from ._utils import set_world_silent
from . import io
from .Topology import Topology
from .DataSetList import DataSetList

_COUNTER = 0

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

    #if not dslist.is_empty():
    #    raise ValueError("do not support non-empty DataSetList")
    if top is None: 
        try: 
           _top = traj.top 
        except: 
            # list, tuple of traj objects 
            _top = traj[0].top 
    elif isinstance(top, string_types):
        _top = Topology(top)
    else:
        _top = top
    if traj is None: 
        raise ValueError("must have trajectory object") 
    elif isinstance(traj, string_types):
        try:
            traj = io.load(traj, _top)
        except:
            raise ValueError("can not load %s" % traj)
    if isinstance(action, string_types): 
        # convert to action 
        act = adict[action] 
    else: 
        act = action 
    return act(command, traj, _top, dslist=dslist, quick_get=quick_get, **kwd)
