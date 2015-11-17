'''sandbox for Julialang and other stuff.
Codes in this module might not be tested and they might be not correct, at all.
'''
import numpy as np
from functools import wraps


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
    indices = traj.top.select(mask)

    for frame in traj:
        diff = np.asarray(to_point) - np.asarray(anchor)
        frame.xyz[indices] += diff
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

    command = ' '.join((filename, mode))
    fi = traj if frame_indices is None else traj.iterframe(
        frame_indices=frame_indices)

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
        command = ' '.join((filename, mode))
        self._outtraj.read_input(command, top=top)
        self._outtraj.process(top)

    def write_frame(self, frame):
        self._outtraj.do_action(frame)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

# for testing
from pytraj._get_common_objects import (_get_data_from_dtype, _get_topology,
                                        _get_reference_from_traj,
                                        _get_fiterator)
from pytraj.actions import CpptrajActions
from pytraj.datasets import CpptrajDatasetList
from pytraj.compat import string_types
from pytraj.utils.convert import array_to_cpptraj_atommask


def _dispatch_traj_ref_top_frame_indices(f):
    @wraps(f)
    def inner(*args, **kwd):
        args = list(args)
        traj = kwd.get('traj', args[0])
        frame_indices = kwd.get('frame_indices')
        ref = kwd.get('ref', None)
        top = kwd.get('top', None)

        if 'mask' in kwd.keys():
            mask = kwd.get('mask')
        else:
            mask = args[1]

        # overwrite
        kwd['top'] = _get_topology(traj, top)
        if ref is not None:
            kwd['ref'] = _get_reference_from_traj(traj, ref)
        if 'traj' in kwd.keys():
            kwd['traj'] = _get_fiterator(traj, frame_indices)
        else:
            args[0] = _get_fiterator(traj, frame_indices)
        if not isinstance(mask, string_types):
            mask = array_to_cpptraj_atommask(mask)
        if 'mask' in kwd.keys():
            kwd['mask'] = mask
        else:
            args[1] = mask
        return f(*args, **kwd)

    return inner


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


def calc_linear_interaction_energy(traj=None,
                                   mask="",
                                   top=None,
                                   dtype='dataset',
                                   frame_indices=None,
                                   *args,
                                   **kwd):
    command = mask
    act = CpptrajActions.Action_LIE()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)

# alias
calc_LIE = calc_linear_interaction_energy

def _dbscan(traj=None,
            mask='*',
            minpoints=None,
            epsilon=0.,
            sievetoframe=False,
            random_sieveseed=1,
            kdist=None,
            kfile=None,
            sieve=1,
            metric='rms',
            top=None,
            options=''):
    '''perform clustering and return cluster index for each frame

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str, default: * (all atoms)
    n_clusters: int, default: 10
    random_point : bool, default: True
    maxit : int, default: 100
        max iterations
    metric : str, {'rms', 'dme'}
        distance metric
    top : Topology, optional, default: None
        only need to provide this Topology if ``traj`` does not have one
    options : option to save data to files.


    Returns
    -------
    1D numpy array of frame indices
    '''

    # don't need to _get_topology
    _top = _get_topology(traj, top)
    _clusters = 'dbscan minpoints ' + str(minpoints)
    _mask = mask
    _epsilon = 'epsilon ' + str(epsilon)
    _sievetoframe = 'sievetoframe' if sievetoframe else ''
    _kdist = 'kdist' + str(kdist) if kdist is not None else ''
    _kfile = kfile if kfile is not None else ''
    _sieve = 'sieve ' + str(sieve)
    _metric = metric
    _random_sieveseed = 'random ' + str(random_sieveseed)
    _output = options
    command = ' '.join((_clusters, _epsilon, _sievetoframe, _kdist, _sieve,
                        _kfile, _metric, _random_sieveseed, _mask, _output))
    return _cluster(traj, command, top=_top, dtype='ndarray')

