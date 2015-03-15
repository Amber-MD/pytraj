from pytraj.externals.six import string_types

def calculate(action=None, command=None, traj=None, top=None, **kwd): 
    # TODO : should write universal help's method 
    """ 
    quick way to get data 
    Parameters: 
    action : Action object or str, default=None 
    command : str, default=None 
    traj : Trajectory object (FrameArray, TrajReadOnly, ...) or list, tuple of traj object 
    top : topology 
 
    Use `calculate(ahelp=True)` or `calculate(ahelp='action name')` for help 
 
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
            adict[kwd['key'].lower()].help() 
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
        return act(command, traj, top, quick_get=True)
