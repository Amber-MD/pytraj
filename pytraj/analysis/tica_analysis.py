"""
Time-Independent Component Analysis (TICA) implementation
"""
import numpy as np
from ..utils.get_common_objects import (
    get_topology,
    get_data_from_dtype,
    get_fiterator,
    super_dispatch,
)
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from ..analysis.c_analysis import c_analysis


__all__ = ['tica']


@super_dispatch()
def tica(traj=None,
         mask='@CA',
         lag=1,
         n_components=None,
         evector_scale='none',
         dtype='dataset',
         top=None,
         frame_indices=None,
         debug=0):
    """Perform Time-Independent Component Analysis (TICA)

    TICA finds the slowest-relaxing collective coordinates in molecular dynamics
    simulations by analyzing time-lagged correlations. This method is particularly
    useful for identifying slow conformational changes and transition pathways.

    Parameters
    ----------
    traj : Trajectory-like
        Input trajectory for TICA analysis
    mask : str, default '@CA'
        Atom selection mask for coordinate extraction
    lag : int, default 1
        Time lag (in frames) for correlation analysis
    n_components : int, optional
        Number of TICA components to compute. If None, compute all possible components
    evector_scale : str, default 'none'
        Eigenvector scaling method:
        - 'none': No scaling
        - 'kinetic': Scale by eigenvalues (kinetic map)
        - 'commute': Scale by regularized time scales (commute map)
    dtype : str, default 'dataset'
        Return data type
    top : Topology, optional
        Topology object
    frame_indices : array-like, optional
        Specific frame indices to analyze
    debug : int, default 0
        Debug level for detailed output

    Returns
    -------
    DatasetList
        TICA results including:
        - Eigenvalues (time scales)
        - Eigenvectors (TICA components)
        - Cumulative variance
        - Projected coordinates (if available)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load('trajectory.nc', 'topology.prmtop')
    >>>
    >>> # Basic TICA analysis on CA atoms with lag=10
    >>> tica_data = pt.tica(traj, mask='@CA', lag=10)
    >>>
    >>> # Get first few components with kinetic scaling
    >>> tica_data = pt.tica(traj, mask='@CA', lag=5, n_components=10,
    ...                     evector_scale='kinetic')
    >>>
    >>> # TICA on backbone heavy atoms
    >>> tica_data = pt.tica(traj, mask='@C,CA,N,O', lag=20)

    Notes
    -----
    - TICA is most effective for identifying slow conformational transitions
    - Larger lag times capture slower dynamics but reduce statistical accuracy
    - The number of meaningful components is often much smaller than the input dimensions
    - For large systems, consider using a subset of atoms (e.g., CA atoms) for efficiency

    References
    ----------
    .. [1] Schwantes, C. R. & Pande, V. S. Improvements in Markov State Model
           Construction Reveal Many Non-Native Contacts in the Folding of NTL9.
           J. Chem. Theory Comput. 9, 2000-2009 (2013).
    .. [2] Pérez-Hernández, G., Paul, F., Giorgino, T., De Fabritiis, G. & Noé, F.
           Identification of slow molecular order parameters for Markov state models.
           J. Chem. Phys. 139, 015102 (2013).
    """

    # Validate parameters
    if lag < 1:
        raise ValueError("lag must be >= 1")

    if evector_scale not in ['none', 'kinetic', 'commute']:
        raise ValueError("evector_scale must be 'none', 'kinetic', or 'commute'")

    # Get trajectory and topology
    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)

    # Build command string for cpptraj TICA analysis
    command_parts = [f"lag {lag}"]

    # Add mask for coordinate selection
    if mask:
        command_parts.append(f"crdset {mask}")

    # Add number of components if specified
    if n_components is not None:
        command_parts.append(f"modes {n_components}")

    # Add eigenvector scaling
    if evector_scale == 'kinetic':
        command_parts.append("kineticmap")
    elif evector_scale == 'commute':
        command_parts.append("commutemap")

    # Add debug level
    if debug > 0:
        command_parts.append(f"debug {debug}")

    # Output datasets
    command_parts.extend([
        "evalout tica_eigenvals.dat",
        "evecout tica_eigenvecs.dat",
        "cumvarout tica_cumvar.dat"
    ])

    command = " ".join(command_parts)

    # Set up TICA analysis
    c_dslist = CpptrajDatasetList()
    analysis = c_analysis.Analysis_TICA()

    try:
        # Initialize and setup analysis
        analysis.read_input(command, dslist=c_dslist)
        analysis.setup(traj.n_frames, topology)

        # Run analysis
        analysis.analyze()

    except Exception as e:
        raise RuntimeError(f"TICA analysis failed: {str(e)}")

    return get_data_from_dtype(c_dslist, dtype)


def tica_transform(traj=None, eigenvectors=None, mask='@CA', top=None, frame_indices=None):
    """Transform trajectory coordinates into TICA space using pre-computed eigenvectors

    Parameters
    ----------
    traj : Trajectory-like
        Input trajectory to transform
    eigenvectors : array-like
        Pre-computed TICA eigenvectors from previous tica() analysis
    mask : str, default '@CA'
        Atom selection mask (must match the one used for computing eigenvectors)
    top : Topology, optional
        Topology object
    frame_indices : array-like, optional
        Specific frame indices to transform

    Returns
    -------
    ndarray
        Transformed coordinates in TICA space

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>>
    >>> # Compute TICA on training trajectory
    >>> train_traj = pt.load('train.nc', 'topology.prmtop')
    >>> tica_data = pt.tica(train_traj, mask='@CA', lag=10)
    >>>
    >>> # Extract eigenvectors (this is conceptual - actual extraction depends on dataset format)
    >>> eigenvectors = tica_data[1].values  # Assuming eigenvectors are in second dataset
    >>>
    >>> # Transform test trajectory
    >>> test_traj = pt.load('test.nc', 'topology.prmtop')
    >>> tica_coords = pt.tica_transform(test_traj, eigenvectors, mask='@CA')
    """

    if eigenvectors is None:
        raise ValueError("eigenvectors must be provided")

    # Get trajectory and topology
    traj = get_fiterator(traj, frame_indices)
    topology = get_topology(traj, top)

    # Get coordinates for the specified mask
    coords = []
    for frame in traj:
        indices = topology.select(mask)
        frame_coords = []
        for idx in indices:
            frame_coords.extend([frame[idx*3], frame[idx*3+1], frame[idx*3+2]])
        coords.append(frame_coords)

    coords = np.array(coords)

    # Center coordinates (subtract mean)
    coords_centered = coords - np.mean(coords, axis=0)

    # Transform using eigenvectors
    if hasattr(eigenvectors, 'values'):
        # If it's a Dataset object
        evecs = eigenvectors.values
    else:
        evecs = np.array(eigenvectors)

    # Project coordinates onto TICA space
    tica_coords = np.dot(coords_centered, evecs)

    return tica_coords


def tica_msm_features(traj=None, mask='@CA', lag=10, n_components=10, top=None):
    """Compute TICA features optimized for Markov State Model construction

    This is a convenience function that combines TICA analysis with best practices
    for MSM construction, including appropriate lag time selection and component filtering.

    Parameters
    ----------
    traj : Trajectory-like
        Input trajectory
    mask : str, default '@CA'
        Atom selection for feature extraction
    lag : int, default 10
        TICA lag time (should be comparable to MSM lag time)
    n_components : int, default 10
        Number of TICA components to retain
    top : Topology, optional
        Topology object

    Returns
    -------
    dict
        Dictionary containing:
        - 'features': TICA-transformed coordinates
        - 'eigenvalues': TICA eigenvalues (time scales)
        - 'cumulative_variance': Cumulative variance explained

    Examples
    --------
    >>> import pytraj as pt
    >>>
    >>> traj = pt.load('md.nc', 'system.prmtop')
    >>> msm_data = pt.tica_msm_features(traj, lag=20, n_components=5)
    >>>
    >>> # Use features for MSM construction
    >>> features = msm_data['features']
    >>> print(f"Feature shape: {features.shape}")
    >>> print(f"Time scales: {msm_data['eigenvalues'][:5]}")
    """

    # Perform TICA analysis
    tica_data = tica(traj, mask=mask, lag=lag, n_components=n_components,
                     evector_scale='kinetic', dtype='dataset', top=top)

    # Extract and organize results
    result = {}

    # Find datasets by name/type
    for i, dataset in enumerate(tica_data):
        key = dataset.key.lower()
        if 'eigenval' in key or 'eval' in key:
            result['eigenvalues'] = dataset.values
        elif 'cumvar' in key or 'variance' in key:
            result['cumulative_variance'] = dataset.values
        elif 'eigenvec' in key or 'evec' in key:
            # Transform trajectory using eigenvectors
            result['features'] = tica_transform(traj, dataset, mask=mask, top=top)

    return result