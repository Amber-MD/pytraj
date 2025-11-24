"""
Correlation and time series analysis functions
"""
import numpy as np
from typing import Union

from ..utils.get_common_objects import (
    get_topology, get_fiterator, get_data_from_dtype,
    get_list_of_commands, get_reference, super_dispatch
)
from ..utils import ensure_not_none_or_string
from ..utils.decorators import register_pmap, register_openmp
from ..utils.convert import array_to_cpptraj_atommask
from ..datasets.datasetlist import DatasetList
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from ..trajectory.shared_methods import iterframe_master
from ..analysis.c_action import c_action, do_action
from ..analysis.c_action.actionlist import ActionList
from .base_classes import (
    CommandType, _check_command_type, _create_and_compute_action_list,
    DatasetType
)


@register_pmap
@super_dispatch()
def atomiccorr(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
               cut=1.0, byres=False):
    """compute atomic position correlation

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for correlation calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    cut : float, default 1.0
        distance cutoff for correlation
    byres : bool, default False
        if True, calculate by residue

    Returns
    -------
    2D ndarray
        correlation matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate atomic correlation for CA atoms
    >>> data = pt.atomiccorr(traj, '@CA')
    >>> # Calculate by residue
    >>> data = pt.atomiccorr(traj, '@CA', byres=True)
    >>> # Use custom cutoff
    >>> data = pt.atomiccorr(traj, '@CA', cut=2.0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['atomiccorr', mask, f'cut {cut}']
    if byres:
        command_parts.append('byres')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_AtomicCorr(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def timecorr(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
             tcorr=10000, norm=False, drct=False):
    """compute time correlation function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for time correlation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    tcorr : int, default 10000
        correlation time
    norm : bool, default False
        if True, normalize the correlation function
    drct : bool, default False
        if True, calculate direct correlation

    Returns
    -------
    1D ndarray
        time correlation function

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate time correlation for coordinates
    >>> data = pt.timecorr(traj, '@CA')
    >>> # Calculate with normalization
    >>> data = pt.timecorr(traj, '@CA', norm=True)
    >>> # Use shorter correlation time
    >>> data = pt.timecorr(traj, '@CA', tcorr=5000)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['timecorr', mask, f'tcorr {tcorr}']
    if norm:
        command_parts.append('norm')
    if drct:
        command_parts.append('drct')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_TimeCorr(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def acorr(traj=None, mask='', frame_indices=None, top=None, dtype='ndarray',
          maxlag=None, norm=False):
    """compute autocorrelation function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for autocorrelation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    maxlag : int, optional
        maximum lag time
    norm : bool, default False
        if True, normalize the autocorrelation

    Returns
    -------
    1D ndarray
        autocorrelation function

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate autocorrelation for distance
    >>> distances = pt.distance(traj, '@1 @10')
    >>> data = pt.acorr(distances)
    >>> # Use custom max lag
    >>> data = pt.acorr(distances, maxlag=100)
    >>> # Normalize
    >>> data = pt.acorr(distances, norm=True)
    """
    ensure_not_none_or_string(traj)

    command_parts = ['autocorr', mask]
    if maxlag is not None:
        command_parts.extend(['maxlag', str(maxlag)])
    if norm:
        command_parts.append('norm')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_AutoCorr,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def xcorr(traj=None, mask1='', mask2='', frame_indices=None, top=None, dtype='ndarray',
          maxlag=None, norm=False):
    """compute cross-correlation function

    Parameters
    ----------
    traj : Trajectory-like
    mask1 : str
        first atom mask
    mask2 : str
        second atom mask
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    maxlag : int, optional
        maximum lag time
    norm : bool, default False
        if True, normalize the cross-correlation

    Returns
    -------
    1D ndarray
        cross-correlation function

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate cross-correlation between two distances
    >>> dist1 = pt.distance(traj, '@1 @5')
    >>> dist2 = pt.distance(traj, '@10 @15')
    >>> data = pt.xcorr(dist1, dist2)
    >>> # Use custom parameters
    >>> data = pt.xcorr(dist1, dist2, maxlag=50, norm=True)
    """
    ensure_not_none_or_string(traj)

    command_parts = ['crosscorr', mask1, mask2]
    if maxlag is not None:
        command_parts.extend(['maxlag', str(maxlag)])
    if norm:
        command_parts.append('norm')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_CrossCorr,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def velocity_autocorrelation(traj=None, mask='', frame_indices=None, top=None,
                           dtype='ndarray', maxlag=None, norm=True,
                           tstep=1.0, usevelocity=False):
    """compute velocity autocorrelation function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for velocity calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    maxlag : int, optional
        maximum lag time
    norm : bool, default True
        if True, normalize the autocorrelation
    tstep : float, default 1.0
        time step between frames
    usevelocity : bool, default False
        if True, use existing velocity data

    Returns
    -------
    1D ndarray
        velocity autocorrelation function

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Calculate velocity autocorrelation
    >>> data = pt.velocity_autocorrelation(traj, '@CA')
    >>> # Use custom time step
    >>> data = pt.velocity_autocorrelation(traj, '@CA', tstep=0.5)
    >>> # Use custom max lag
    >>> data = pt.velocity_autocorrelation(traj, '@CA', maxlag=200)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['velocityautocorr', mask, f'tstep {tstep}']
    if maxlag is not None:
        command_parts.extend(['maxlag', str(maxlag)])
    if norm:
        command_parts.append('norm')
    if usevelocity:
        command_parts.append('usevelocity')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_VelocityAutoCorr(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


# Create alias
velocityautocorr = velocity_autocorrelation


def auto_correlation_function(data, maxlag=None, norm=False):
    """compute autocorrelation function for numpy array

    Parameters
    ----------
    data : array-like
        input data series
    maxlag : int, optional
        maximum lag time (default: len(data)//4)
    norm : bool, default False
        if True, normalize the autocorrelation

    Returns
    -------
    1D ndarray
        autocorrelation function

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> # Generate some test data
    >>> data = np.random.randn(1000)
    >>> acorr = pt.auto_correlation_function(data)
    >>> # Use custom parameters
    >>> acorr = pt.auto_correlation_function(data, maxlag=100, norm=True)
    """
    data = np.asarray(data)
    n = len(data)

    if maxlag is None:
        maxlag = n // 4

    maxlag = min(maxlag, n - 1)

    # Remove mean
    data_centered = data - np.mean(data)

    # Compute autocorrelation using FFT
    # Pad with zeros to twice the length
    padded = np.zeros(2 * n)
    padded[:n] = data_centered

    # FFT, multiply by conjugate, then inverse FFT
    fft_data = np.fft.fft(padded)
    autocorr_full = np.fft.ifft(fft_data * np.conj(fft_data)).real

    # Take only the first half and normalize by decreasing sample size
    autocorr = autocorr_full[:maxlag + 1]
    for i in range(len(autocorr)):
        autocorr[i] /= (n - i)

    if norm:
        autocorr /= autocorr[0]

    return autocorr


def cross_correlation_function(data1, data2, maxlag=None, norm=False):
    """compute cross-correlation function for two numpy arrays

    Parameters
    ----------
    data1 : array-like
        first data series
    data2 : array-like
        second data series
    maxlag : int, optional
        maximum lag time (default: min(len(data1), len(data2))//4)
    norm : bool, default False
        if True, normalize the cross-correlation

    Returns
    -------
    1D ndarray
        cross-correlation function

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> # Generate test data
    >>> data1 = np.random.randn(1000)
    >>> data2 = np.random.randn(1000)
    >>> xcorr = pt.cross_correlation_function(data1, data2)
    >>> # Use custom parameters
    >>> xcorr = pt.cross_correlation_function(data1, data2, maxlag=100, norm=True)
    """
    data1 = np.asarray(data1)
    data2 = np.asarray(data2)

    n1, n2 = len(data1), len(data2)
    n = min(n1, n2)

    if maxlag is None:
        maxlag = n // 4

    maxlag = min(maxlag, n - 1)

    # Remove means
    data1_centered = data1[:n] - np.mean(data1[:n])
    data2_centered = data2[:n] - np.mean(data2[:n])

    # Compute cross-correlation using FFT
    # Pad with zeros
    padded1 = np.zeros(2 * n)
    padded2 = np.zeros(2 * n)
    padded1[:n] = data1_centered
    padded2[:n] = data2_centered

    # FFT, multiply data1 by conjugate of data2, then inverse FFT
    fft1 = np.fft.fft(padded1)
    fft2 = np.fft.fft(padded2)
    xcorr_full = np.fft.ifft(fft1 * np.conj(fft2)).real

    # Take correlation at positive lags and normalize
    xcorr = xcorr_full[:maxlag + 1]
    for i in range(len(xcorr)):
        xcorr[i] /= (n - i)

    if norm:
        # Normalize by the geometric mean of autocorrelations at zero lag
        auto1_0 = np.sum(data1_centered**2) / n
        auto2_0 = np.sum(data2_centered**2) / n
        xcorr /= np.sqrt(auto1_0 * auto2_0)

    return xcorr


def coherent_neutron_scattering(traj=None, mask='', frame_indices=None, top=None,
                               dtype='ndarray', qmax=50.0, nq=50):
    """compute coherent neutron scattering function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for scattering calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    qmax : float, default 50.0
        maximum q value
    nq : int, default 50
        number of q points

    Returns
    -------
    DatasetList
        scattering data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.coherent_neutron_scattering(traj, ':1-13')
    """
    ensure_not_none_or_string(traj)

    command = f'cns {mask} qmax {qmax} nq {nq}'
    dslist, _ = do_action(traj, command, c_action.Action_CNS,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def incoherent_neutron_scattering(traj=None, mask='', frame_indices=None, top=None,
                                 dtype='ndarray', qmax=50.0, nq=50):
    """compute incoherent neutron scattering function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for scattering calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    qmax : float, default 50.0
        maximum q value
    nq : int, default 50
        number of q points

    Returns
    -------
    DatasetList
        scattering data

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.incoherent_neutron_scattering(traj, '@H*')
    """
    ensure_not_none_or_string(traj)

    command = f'ined {mask} qmax {qmax} nq {nq}'
    dslist, _ = do_action(traj, command, c_action.Action_INED,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'atomiccorr', 'acorr', 'xcorr', 'timecorr', 'velocity_autocorrelation',
    'velocityautocorr', 'auto_correlation_function', 'cross_correlation_function',
    'coherent_neutron_scattering', 'incoherent_neutron_scattering'
]