"""
Statistical analysis functions
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
def pca(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
        altmask='', ref=None, bymass=False):
    """Principal Component Analysis

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for PCA calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    altmask : str, optional
        alternative mask for different coordinate set
    ref : Frame or int, optional
        reference structure for alignment
    bymass : bool, default False
        if True, mass-weight the coordinates

    Returns
    -------
    DatasetList
        PCA results including eigenvalues, eigenvectors, and projections

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # PCA on CA atoms
    >>> pca_data = pt.pca(traj, '@CA')
    >>> # Mass-weighted PCA
    >>> pca_data = pt.pca(traj, '@CA', bymass=True)
    >>> # PCA with reference structure
    >>> pca_data = pt.pca(traj, '@CA', ref=0)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    command_parts = ['matrix', 'covar', mask]
    if altmask:
        command_parts.extend(['altmask', altmask])
    if bymass:
        command_parts.append('bymass')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Matrix(), command, top, dslist=dslist)
    action_list.compute(traj)

    # Diagonalize the covariance matrix
    diagmatrix_command = 'diagmatrix out evecs.dat vecs 10'
    action_list.add(c_action.Action_DiagMatrix(), diagmatrix_command, top, dslist=dslist)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def projection(traj=None, mask='', modes=None, frame_indices=None, top=None,
               dtype='ndarray', bymass=False):
    """project trajectory onto principal components or normal modes

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for projection
    modes : array-like or str
        eigenvectors/modes for projection
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'ndarray'
    bymass : bool, default False
        if True, mass-weight the projection

    Returns
    -------
    2D ndarray
        projection data, shape (n_frames, n_modes)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # First do PCA to get modes
    >>> pca_data = pt.pca(traj, '@CA')
    >>> # Project trajectory onto first 3 modes
    >>> proj = pt.projection(traj, '@CA', modes='evecs.dat[:3]')
    >>> # Mass-weighted projection
    >>> proj = pt.projection(traj, '@CA', modes='evecs.dat[:3]', bymass=True)
    """
    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top = get_topology(traj, top)

    if modes is None:
        raise ValueError("modes parameter is required for projection")

    command_parts = ['projection', mask]
    if isinstance(modes, str):
        command_parts.extend(['modes', modes])
    if bymass:
        command_parts.append('bymass')

    command = ' '.join(command_parts)

    dslist = CpptrajDatasetList()
    action_list = ActionList()
    action_list.add(c_action.Action_Projection(), command, top, dslist=dslist)
    action_list.compute(traj)

    return get_data_from_dtype(dslist, dtype)


@register_pmap
@super_dispatch()
def wavelet(traj=None, mask='', frame_indices=None, top=None, dtype='dataset',
            nb=0, s0=1.0, ds=0.25, correction=1.0):
    """wavelet analysis of time series

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for wavelet analysis
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    nb : int, default 0
        number of points to add as buffer
    s0 : float, default 1.0
        smallest scale
    ds : float, default 0.25
        scale spacing
    correction : float, default 1.0
        correction factor

    Returns
    -------
    DatasetList
        wavelet analysis results

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Wavelet analysis of RMSD time series
    >>> rmsd_data = pt.rmsd(traj, '@CA')
    >>> wavelet_data = pt.wavelet(rmsd_data)
    >>> # Custom wavelet parameters
    >>> wavelet_data = pt.wavelet(rmsd_data, s0=0.5, ds=0.1)
    """
    ensure_not_none_or_string(traj)

    command_parts = ['wavelet', mask, f'nb {nb}', f's0 {s0}',
                    f'ds {ds}', f'correction {correction}']

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Wavelet,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def covariance_matrix(traj=None, mask='', frame_indices=None, top=None,
                     dtype='dataset', bymass=False):
    """compute covariance matrix

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for covariance calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'
    bymass : bool, default False
        if True, mass-weight the covariance

    Returns
    -------
    DatasetList
        covariance matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Covariance matrix for CA atoms
    >>> cov_matrix = pt.covariance_matrix(traj, '@CA')
    >>> # Mass-weighted covariance
    >>> cov_matrix = pt.covariance_matrix(traj, '@CA', bymass=True)
    """
    ensure_not_none_or_string(traj)

    command_parts = ['matrix', 'covar', mask]
    if bymass:
        command_parts.append('bymass')

    command = ' '.join(command_parts)

    dslist, _ = do_action(traj, command, c_action.Action_Matrix,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def correlation_matrix(traj=None, mask='', frame_indices=None, top=None,
                      dtype='dataset'):
    """compute correlation matrix

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for correlation calculation
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList
        correlation matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Correlation matrix for CA atoms
    >>> corr_matrix = pt.correlation_matrix(traj, '@CA')
    """
    ensure_not_none_or_string(traj)

    command = f'matrix correl {mask}'

    dslist, _ = do_action(traj, command, c_action.Action_Matrix,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


def distance_covariance(traj=None, mask='', frame_indices=None, top=None,
                       dtype='dataset'):
    """compute distance covariance matrix

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
        atom mask for distance covariance
    frame_indices : array-like, optional, default None
    top : Topology, optional
    dtype : str, default 'dataset'

    Returns
    -------
    DatasetList
        distance covariance matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # Distance covariance for CA atoms
    >>> dist_cov = pt.distance_covariance(traj, '@CA')
    """
    ensure_not_none_or_string(traj)

    command = f'matrix distcovar {mask}'

    dslist, _ = do_action(traj, command, c_action.Action_Matrix,
                         frame_indices=frame_indices, top=top)
    return get_data_from_dtype(dslist, dtype)


# Export all functions
__all__ = [
    'pca', 'projection', 'wavelet', 'covariance_matrix',
    'correlation_matrix', 'distance_covariance'
]