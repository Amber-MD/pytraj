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

class Trajout:
    # give wrong n_atoms for last frame. Why?
    '''
    Examples
    --------

    >>> with Trajout('test.pdb', mode='model', top=top) trajout:
    >>>     for idx, frame in enumerate(traj):
    >>>         trajout.write_frame(idx, frame)

    '''
    def __init__(self, filename, mode='', top=None):
        from pytraj.actions.CpptrajActions import Action_Outtraj
        self._outtraj = Action_Outtraj()
        command = ' ' .join((filename, mode))
        self._outtraj.read_input(command, top=top)
        self._outtraj.process(top)

    def write_frame(self, frame):
        self._outtraj.do_action(frame)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass
