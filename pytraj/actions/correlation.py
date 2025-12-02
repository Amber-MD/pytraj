"""
Correlation and time series analysis functions
"""
from .base import *

__all__ = [
    'atomiccorr', 'acorr', 'xcorr', 'timecorr', 'velocity_autocorrelation',
    'velocityautocorr', 'crank', 'wavelet'
]


@super_dispatch()
def atomiccorr(traj,
               mask='',
               mask2='',
               cut=0.0,
               min_spacing=1.0,
               dtype='ndarray',
               byres=False,
               top=None,
               frame_indices=None):
    """compute average correlations between the motion of atoms in mask.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        first atom mask
    mask2 : str, optional
        second atom mask (if not provided, correlations within mask will be computed)
    cut : float, default 0.0
        if greater than 0, only print correlations with absolute value greater than cut
    min_spacing : float, default 1.0
        only calculate correlations for motion vectors spaced min_spacing apart
    dtype : str, default 'ndarray'
        return data type
    byres : bool, default False
        if False, compute atomic motion vector
        if True, calculate motion vectors for entire residues (selected atoms in residues only).
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray or DatasetList
        correlation matrix between atomic motions

    Notes
    -----
    This function computes correlations between atomic motion vectors. The implementation
    has been updated from the original backup version to support two-mask correlations
    and different default parameter values.
    """
    command = mask
    if mask2:
        command += f" {mask2}"
    if cut > 0:
        command += f" cut {cut}"
    if min_spacing != 1.0:
        command += f" min {min_spacing}"
    if byres:
        command += " byres"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_AtomicCorr()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


def timecorr(vec0, vec1, order=2, tstep=1., tcorr=10000., norm=False, dtype='ndarray'):
    """Compute time correlation.

    Parameters
    ----------
    vec0 : 2D array-like, shape=(n_frames, 3)
    vec1 : 2D array-like, shape=(n_frames, 3)
    order : int, default 2
    tstep : float, default 1.
    tcorr : float, default 10000.
    norm : bool, default False
    dtype : str, default 'ndarray'
    """
    runner = AnalysisRunner(c_analysis.Analysis_Timecorr)
    runner.add_dataset(DatasetType.VECTOR, "_vec0", vec0)
    runner.add_dataset(DatasetType.VECTOR, "_vec1", vec1)

    command = f"vec1 _vec0 vec2 _vec1 order {order} tstep {tstep} tcorr {tcorr} {'norm' if norm else ''}"
    runner.run_analysis(command)

    return get_data_from_dtype(runner.datasets[2:], dtype=dtype)


@super_dispatch()
def velocity_autocorrelation(
        traj,
        mask='',
        maxlag=-1,
        tstep=1.0,
        direct=True,
        norm=False,
        usecoords=False,
        dtype='ndarray',
        top=None,
        velocity_arr=None):
    """
    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask
    maxlag : int, default -1
        maximum lag. If -1, using half of total number of frame
        if given, use it.
    tstep : float, default 1.0
    direct : bool, default True
        if True, use direct method
        else, use FFT to compute autocorrelation function
    norm : bool, default False
        if True, normalize autocorrelation function to 1.0
    usecoords : bool, default False
        if True, use velocity info in Frame
    dtype : str, default 'ndarray'
        return data type
    top : None or Topology, default None, optional
    velocity_arr : None or 3D like array, default None
        only use `velocity_arr` if usecoords is True

    Notes
    -----
    If you create Trajectory by `pytraj.load` method, there is no velocity information.
    So if you want to use `usecoords=True`, you need to provide 3D-array velocity_arr
    """
    from pytraj import Frame

    if velocity_arr is not None:
        velocity_arr = np.asarray(velocity_arr)

        if len(velocity_arr.shape) != 3:
            raise ValueError(
                'provided velocity_arr must be 3D array-like, shape=(n_frames, n_atoms, 3)'
            )

    velocity_autocorrelation_action = c_action.Action_VelocityAutoCorr()
    action_datasets = CpptrajDatasetList()

    command = f"maxlag {maxlag} tstep {tstep} {'direct' if direct else ''} {'norm' if norm else ''} {'usecoords' if usecoords else ''}"
    crdinfo = dict(has_velocity=True)

    velocity_autocorrelation_action.read_input(command, top, dslist=action_datasets)
    velocity_autocorrelation_action.setup(top, crdinfo=crdinfo)

    frame_template = Frame()

    if usecoords and velocity_arr is not None:
        frame_template._allocate_force_and_velocity(
            top, crdinfo=dict(has_velocity=True))
        use_template = True
    else:
        use_template = False

    for idx, frame in enumerate(traj):
        if not use_template:
            if usecoords and not frame.has_velocity():
                raise ValueError(
                    "Frame must have velocity if specify 'usecoords'")
            velocity_autocorrelation_action.compute(frame)
        else:
            vel = velocity_arr[idx]
            frame_template.xyz[:] = frame.xyz[:]
            frame_template.velocity[:] = vel
            velocity_autocorrelation_action.compute(frame_template)

    velocity_autocorrelation_action.post_process()

    return get_data_from_dtype(action_datasets, dtype=dtype)


# create an alias
velocityautocorr = velocity_autocorrelation


def crank(data0, data1, mode='distance', dtype='ndarray'):
    """compute Crank-Nicholson reaction coordinate

    Parameters
    ----------
    data0 : array-like
        first dataset
    data1 : array-like
        second dataset
    mode : str, default 'distance'
        calculation mode
    dtype : str, default 'ndarray'
        return data type

    Returns
    -------
    out : ndarray
    """
    data0 = np.asarray(data0, dtype='f8')
    data1 = np.asarray(data1, dtype='f8')

    command = f"crank {mode}"

    c_dslist = CpptrajDatasetList()

    # add datasets
    dataset0 = c_dslist.add('double', 'data0')
    dataset0.data = data0
    dataset1 = c_dslist.add('double', 'data1')
    dataset1.data = data1

    # run analysis
    c_analysis.Analysis_CrankShaft(command, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype=dtype)


def acorr(data, dtype='ndarray', option=''):
    """compute autocorrelation

    Parameters
    ----------
    data : 1d array-like
    dtype: return type, default 'ndarray'
    covar : bool, default True
    option : str
        more cpptraj options

    Notes
    -----
    Same as `autocorr` in cpptraj
    """
    runner = AnalysisRunner(c_analysis.Analysis_AutoCorr)
    runner.add_dataset(DatasetType.DOUBLE, "d0", np.asarray(data))

    command = "d0 out _tmp.out"
    runner.run_analysis(command)

    return get_data_from_dtype(runner.datasets[1:], dtype=dtype)


def xcorr(data0, data1, dtype='ndarray'):
    """compute cross correlation between two datasets

    Parameters
    ----------
    data0 and data1: 1D-array like
    dtype : return datatype, default 'ndarray'


    Notes
    -----
    Same as `corr` in cpptraj
    """
    runner = AnalysisRunner(c_analysis.Analysis_Corr)
    runner.add_dataset(DatasetType.DOUBLE, "d0", np.asarray(data0))
    runner.add_dataset(DatasetType.DOUBLE, "d1", np.asarray(data1))

    command = "d0 d1 out _tmp.out"
    runner.run_analysis(command)

    return get_data_from_dtype(runner.datasets[2:3], dtype=dtype)


def wavelet(traj, command):
    """wavelet analysis

    Parameters
    ----------
    traj : Trajectory-like
    command : str, cpptraj command

    Returns
    -------
    out : dict

    Notes
    -----
    - This method is not well-supported in pytraj. It means that
    you need to type cpptraj command. Please check cpptraj manual for further
    info if you really want to use it.

    - Currently pytraj will create a new copy of Trajectory for cpptraj in memory,
    so this method is only good for small trajectory that fit to your RAM.

    version added: 1.0.6

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_dpdp()
    >>> c0 = 'nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 :1-22'
    >>> c1 = 'cluster minpoints 66 epsilon 10.0'
    >>> command = ' '.join((c0, c1))
    >>> wavelet_dict = pt.wavelet(traj, command)
    """
    runner = AnalysisRunner(c_analysis.Analysis_Wavelet)
    runner.add_dataset(DatasetType.COORDS, "_DEFAULTCRD_", traj)
    runner.run_analysis(command)
    runner.datasets.remove_set(runner.datasets["_DEFAULTCRD_"])
    return get_data_from_dtype(runner.datasets, dtype='dict')