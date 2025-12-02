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
    """compute atomic correlation

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        first mask
    mask2 : str, optional
        second mask
    cut : float, default 0.0
        cutoff distance
    min_spacing : float, default 1.0
        minimum spacing between atoms
    dtype : str, default 'ndarray'
        return data type
    byres : bool, default False
        if True, compute per residue
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray or DatasetList
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
        mask='*',
        tstep=1.0,
        tcorr=10000.0,
        order=2,
        norm=False,
        direct=False,
        dtype='ndarray',
        top=None,
        frame_indices=None):
    """compute velocity autocorrelation <v(0).v(t)>

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*'
        atom mask
    tstep : float, default 1.0
        time step between frames (ps)
    tcorr : float, default 10000.0
        correlation time (ps)
    order : int, default 2
        1: |v(0)|*|v(t)|*cos(theta)  2: v(0).v(t)
    norm : bool, default False
        if True, normalize correlation function
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_time_points,)

    Notes
    -----
    need velocity information
    """
    command = mask + f" tstep {tstep} tcorr {tcorr} order {order}"
    if norm:
        command += " norm"
    if direct:
        command += " direct"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_VelocityAutoCorr()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


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