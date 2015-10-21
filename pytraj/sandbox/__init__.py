'''sandbox for Julialang and other stuff.
Codes in this module might not be tested and they might be not correct, at all.
'''
import numpy as np


def take(traj, indices):
    return traj[indices]


def itake(traj, indices):
    return traj.iterframe(frame_indices=indices)


def get_top(traj):
    return traj.top


def set_top(traj, top):
    traj.top = top

def translate(traj, mask, anchor, to_point=[0.0, 0.0, 0.0]):
    '''translate a group of atoms in mask by moving anchor mask to given point

    Returns
    -------
    updated traj
    '''
    from pytraj import center_of_mass

    indices = traj.top.select(mask)

    for frame in traj:
        center = center_of_mass(frame, mask=anchor, top=traj.top)[0]
        diff = np.asarray(to_point) - center
        frame.xyz += diff
    return traj


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

# for testing
from pytraj._get_common_objects import (_dispatch_traj_ref_top_frame_indices,
        _get_data_from_dtype)
from pytraj.actions import CpptrajActions
from pytraj.datasets import CpptrajDatasetList 

@_dispatch_traj_ref_top_frame_indices
def _toy_radgyr(traj,
                mask="",
                top=None,
                dtype='ndarray',
                nomax=True,
                frame_indices=None):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.mindist(traj, '@CA @H')
    '''
    act = CpptrajActions.Action_Radgyr()
    dslist = CpptrajDatasetList()

    _nomax = 'nomax' if nomax else ''
    mask = ' '.join((mask, _nomax))

    act(mask, traj, top=top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)
