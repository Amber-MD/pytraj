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
    """compute time correlation: <v(i).v(i+t)>

    Parameters
    ----------
    vec0 : array-like, shape=(n_frames, 3)
    vec1 : array-like, shape=(n_frames, 3)
    order : int, default 2
        1 : first order (|v(i)|*|v(i+t)|*cos(theta(i, i+t)))
        2 : second order (v(i) . v(i+t))
    tstep : float, default 1.
        time step between frames
    tcorr : float, default 10000.
        correlation time
    norm : bool, default False
        if True, normalize
    dtype : str, default 'ndarray'
        return data type

    Returns
    -------
    out : ndarray
    """
    vec0 = np.asarray(vec0, dtype='f8')
    vec1 = np.asarray(vec1, dtype='f8')

    if vec0.ndim != 2:
        raise ValueError("vec0 must be 2D array")
    if vec1.ndim != 2:
        raise ValueError("vec1 must be 2D array")
    if vec0.shape != vec1.shape:
        raise ValueError("vec0 and vec1 must have same shape")

    max_frames = int(tcorr / tstep)
    if max_frames >= len(vec0):
        max_frames = len(vec0) - 1

    result = np.zeros(max_frames + 1)

    if order == 1:
        # first order: |v(i)|*|v(i+t)|*cos(theta(i, i+t))
        for t in range(max_frames + 1):
            corr_sum = 0.0
            count = 0
            for i in range(len(vec0) - t):
                v_i = vec0[i]
                v_j = vec1[i + t]
                mag_i = np.linalg.norm(v_i)
                mag_j = np.linalg.norm(v_j)
                if mag_i > 0 and mag_j > 0:
                    cos_theta = np.dot(v_i, v_j) / (mag_i * mag_j)
                    corr_sum += mag_i * mag_j * cos_theta
                    count += 1
            if count > 0:
                result[t] = corr_sum / count
    else:
        # second order: v(i) . v(i+t)
        for t in range(max_frames + 1):
            corr_sum = 0.0
            count = 0
            for i in range(len(vec0) - t):
                corr_sum += np.dot(vec0[i], vec1[i + t])
                count += 1
            if count > 0:
                result[t] = corr_sum / count

    if norm and result[0] != 0:
        result /= result[0]

    return result


@super_dispatch()
def velocity_autocorrelation(
        traj,
        mask='*',
        tstep=1.0,
        tcorr=10000.0,
        order=2,
        norm=False,
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
    data : array-like
        input data
    dtype : str, default 'ndarray'
        return data type
    option : str, optional
        extra options

    Returns
    -------
    out : ndarray
    """
    data = np.asarray(data, dtype='f8')

    c_dslist = CpptrajDatasetList()

    # add dataset
    dataset = c_dslist.add('double', 'data')
    dataset.data = data

    command = f"autocorr data {option}"

    # run analysis
    c_analysis.Analysis_AutoCorr(command, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype=dtype)


def xcorr(data0, data1, dtype='ndarray'):
    """compute cross correlation between two datasets

    Parameters
    ----------
    data0 : array-like
        first dataset
    data1 : array-like
        second dataset
    dtype : str, default 'ndarray'
        return data type

    Returns
    -------
    out : ndarray
    """
    data0 = np.asarray(data0, dtype='f8')
    data1 = np.asarray(data1, dtype='f8')

    c_dslist = CpptrajDatasetList()

    # add datasets
    dataset0 = c_dslist.add('double', 'data0')
    dataset0.data = data0
    dataset1 = c_dslist.add('double', 'data1')
    dataset1.data = data1

    command = "corr data0 data1"

    # run analysis
    c_analysis.Analysis_CrossCorr(command, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype=dtype)


def wavelet(traj, command):
    """perform wavelet analysis

    Parameters
    ----------
    traj : Trajectory-like
    command : str
        cpptraj wavelet command

    Returns
    -------
    out : DatasetList
    """
    c_dslist = CpptrajDatasetList()

    # wavelet is an analysis, not an action
    c_analysis.Analysis_Wavelet(command, dslist=c_dslist)

    return c_dslist