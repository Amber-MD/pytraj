'''sandbox for Julialang and other stuff.
Codes in this module might not be tested and they might be not correct, at all.
'''


def take(traj, indices):
    return traj[indices]


def itake(traj, indices):
    return traj.iterframe(frame_indices=indices)


def get_top(traj):
    return traj.top


def set_top(traj, top):
    traj.top = top

def write_traj(filename, traj=None, mode='', frame_indices=None):
    '''
    Parameters
    ----------

    Notes
    -----
    cpptraj will detect file format based on extension for writting.


    Examples
    --------
    '''
    from pytraj.actions.CpptrajActions import Action_Outtraj

    command = ' ' .join((filename, mode))
    fi = traj if frame_indices is None else traj.iterframe(frame_indices=frame_indices)

    act = Action_Outtraj()
    _top = traj.top
    act(command, fi, top=_top)
