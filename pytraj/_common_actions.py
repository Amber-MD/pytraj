from pytraj.externals.six import string_types
from pytraj._utils import set_world_silent

def calculate(action=None, command=None, traj=None, top=None, **kwd): 
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
    from pytraj import adict 
    if action is None and command is None and traj is None and top is None: 
        if not kwd: 
            # 
            #print (calculate.__doc__) 
            print () 
            print (adict.keys()) 
            print () 
            print ("use calculate(key=action_name) for help") 
        else: 
            set_world_silent(False)
            adict[kwd['key'].lower()].help() 
            set_world_silent(True)
    else: 
        if top is None: 
            try: 
               top = traj.top 
            except: 
                # list, tuple of traj objects 
                top = traj[0].top 
        if traj is None: 
            raise ValueError("must have trajectory object") 
        if isinstance(action, string_types): 
            # convert to action 
            act = adict[action] 
        else: 
            act = action 
        return act(command, traj, top, quick_get=True, **kwd)
