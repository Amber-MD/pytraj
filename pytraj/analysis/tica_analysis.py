"""
Time-Independent Component Analysis (TICA) using cpptraj
"""
import numpy as np
from ..utils.get_common_objects import (
    get_data_from_dtype,
)
from ..datasets.c_datasetlist import DatasetList as CpptrajDatasetList


__all__ = ['tica']


def _is_angular_data(data):
    """Check if data appears to be angular (dihedral) data in degrees"""
    data = np.asarray(data)
    # Heuristic for detecting dihedral angles - adjusted for typical dihedral ranges
    return (
        data.min() >= -180 and data.max() <= 180 and  # Strict dihedral range
        (data.max() - data.min()) > 30 and  # Some angular variation (lowered threshold)
        np.abs(data).max() > 30  # Must have moderate angles (lowered threshold)
    )


def tica(traj=None,
         mask='@CA',
         data=None,
         lag=1,
         n_components=None,
         evector_scale='none',
         dtype='dataset',
         top=None,
         frame_indices=None,
         commute=False,
         cumvarout=None):
    """Perform Time-Independent Component Analysis (TICA) using cpptraj

    TICA finds the slowest-relaxing collective coordinates in molecular dynamics
    simulations by analyzing time-lagged correlations. This method is particularly
    useful for identifying slow conformational changes and transition pathways.

    Parameters
    ----------
    traj : Trajectory-like, optional
        Input trajectory for coordinate-based TICA analysis
    mask : str, default '@CA'
        Atom selection mask for coordinate extraction (used with traj)
    data : list of arrays, optional
        Pre-computed data arrays for dataset-based TICA analysis.
        Alternative to using traj+mask. Each array should be 1D with same length.
    lag : int, default 1
        Time lag (in frames) for correlation analysis
    n_components : int, optional
        Number of TICA components to compute. If None, compute all possible components
    evector_scale : str, default 'none'
        Eigenvector scaling method (passed to cpptraj)
    dtype : str, default 'dataset'
        Return data type
    top : Topology, optional
        Topology object
    frame_indices : array-like, optional
        Specific frame indices to analyze

    Returns
    -------
    CpptrajDatasetList
        TICA results from cpptraj's Analysis_TICA implementation

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load('trajectory.nc', 'topology.prmtop')
    >>>
    >>> # Basic TICA analysis on CA atoms with lag=10
    >>> tica_data = pt.tica(traj, mask='@CA', lag=10)
    >>>
    >>> # Dataset-based TICA with distances
    >>> d1 = pt.distance(traj, ':1@CA :5@CA')
    >>> d2 = pt.distance(traj, ':2@CA :6@CA')
    >>> tica_data = pt.tica(data=[d1, d2], lag=5)
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
    from ..utils.get_common_objects import get_iterator_from_dslist, get_data_from_dtype
    from ..analysis.c_analysis.c_analysis import Analysis_TICA

    # Validate parameters
    if lag < 1:
        raise ValueError("lag must be >= 1")

    # Use cpptraj's Analysis_TICA implementation
    if data is not None and traj is not None:
        raise ValueError("Cannot specify both 'data' and 'traj'. Use one or the other.")
    elif data is not None:
        return _tica_dataset_based_cpptraj(data, lag, n_components, evector_scale, dtype, commute, cumvarout)
    elif traj is not None:
        return _tica_coordinate_based_cpptraj(traj, mask, lag, n_components, evector_scale, dtype, top, frame_indices, commute, cumvarout)
    else:
        raise ValueError("Must specify either 'data' or 'traj' parameter.")


def _tica_coordinate_based_cpptraj(traj, mask, lag, n_components, evector_scale, dtype, top, frame_indices, commute, cumvarout):
    """Coordinate-based TICA using cpptraj Analysis_TICA"""
    from ..actions.base import AnalysisRunner, DatasetType, CommandBuilder
    from ..analysis.c_analysis.c_analysis import Analysis_TICA
    from ..utils.get_common_objects import get_data_from_dtype

    # Use AnalysisRunner pattern for coordinate-based analysis
    runner = AnalysisRunner(Analysis_TICA)

    # Apply frame_indices if specified
    if frame_indices is not None:
        subset_traj = traj[frame_indices]
    else:
        subset_traj = traj

    # Add trajectory as coordinate dataset
    runner.add_dataset(DatasetType.COORDS, "_DEFAULTCRD_", subset_traj)

    # Build TICA command for coordinate-based analysis using CommandBuilder
    command = (CommandBuilder()
               .add("crdset", "_DEFAULTCRD_")
               .add("mask", mask)
               .add("lag", lag)
               .add("evecs", n_components, condition=n_components is not None)
               .add("scale", evector_scale, condition=evector_scale != 'none')
               .add("map commute", condition=commute)
               .add("cumvarout", cumvarout, condition=cumvarout is not None)
               .add("out tica_output.dat", condition=cumvarout is None)
               .build())

    # Run the analysis
    runner.run_analysis(command)

    return _extract_tica_results(runner.datasets, dtype, n_components)


def _tica_dataset_based_cpptraj(data, lag, n_components, evector_scale, dtype, commute, cumvarout):
    """Dataset-based TICA using cpptraj Analysis_TICA"""
    from ..actions.base import CommandBuilder
    from ..analysis.c_analysis.c_analysis import Analysis_TICA
    from ..utils.get_common_objects import get_data_from_dtype

    # Create cpptraj Analysis_TICA object
    act = Analysis_TICA()

    # Create dataset list and add data arrays
    c_dslist = CpptrajDatasetList()

    if isinstance(data, (list, tuple)):
        processed_datasets = []

        # Process each dataset - convert angular data to sin/cos components
        for i, dataset in enumerate(data):
            dataset_array = np.asarray(dataset, dtype='f8')

            # Check if this looks like angular data (degrees, typical dihedral range)
            if _is_angular_data(dataset_array):
                # Convert to radians and create sin/cos components
                rad_data = np.deg2rad(dataset_array)
                sin_component = np.sin(rad_data)
                cos_component = np.cos(rad_data)
                processed_datasets.extend([sin_component, cos_component])
            else:
                # Keep non-angular data as-is
                processed_datasets.append(dataset_array)

        # Add processed datasets to cpptraj dataset list
        for i, proc_data in enumerate(processed_datasets):
            name = f'data_{i}'
            c_dslist.add('double', name)
            c_dslist[-1].data = proc_data

        # Build data command referencing all datasets
        data_refs = ' '.join([f'data data_{i}' for i in range(len(processed_datasets))])
    else:
        raise ValueError("Data must be a list or tuple of arrays for dataset-based TICA")

    # Build TICA command for cpptraj using CommandBuilder
    tica_command = (CommandBuilder()
                    .add(data_refs, condition=True)
                    .add("lag", lag)
                    .add("evecs", n_components, condition=n_components is not None)
                    .add("scale", evector_scale, condition=evector_scale != 'none')
                    .add("map commute", condition=commute)
                    .add("cumvarout", cumvarout, condition=cumvarout is not None)
                    .add("out tica_output.dat", condition=cumvarout is None)
                    .build())

    # Run cpptraj TICA analysis
    act(tica_command, dslist=c_dslist)

    return _extract_tica_results(c_dslist, dtype, n_components)


def _extract_tica_results(c_dslist, dtype, n_components=None):
    """Extract TICA results from cpptraj dataset list and return in expected format"""

    # Create a simple result object that mimics the old TICAResult class
    class TICAResult:
        def __init__(self):
            self.cumvar = None
            self.eigenvalues = None
            self.eigenvectors = None
            self.modes = None  # Add modes attribute

    result = TICAResult()

    # Find and extract the key datasets from cpptraj output
    for dataset in c_dslist:
        if hasattr(dataset, 'key'):
            key = dataset.key

            # Look for cumulative variance data (ends with [cumvar])
            if '[cumvar]' in key:
                if hasattr(dataset, 'values'):
                    result.cumvar = np.asarray(dataset.values).copy()

            # Look for DatasetModes which contains eigenvalues and eigenvectors
            elif 'TICA_' in key and not '[' in key:  # Main TICA dataset without brackets
                if hasattr(dataset, 'values'):
                    try:
                        modes_data = dataset.values
                        # modes_data should be a tuple: (eigenvalues, eigenvectors)
                        if isinstance(modes_data, tuple) and len(modes_data) == 2:
                            result.eigenvalues = np.asarray(modes_data[0])
                            result.eigenvectors = np.asarray(modes_data[1])
                            # modes should have shape (n_features, n_components)
                            result.modes = result.eigenvectors if result.eigenvectors is not None else None

                            # Apply n_components slicing if specified
                            if n_components is not None:
                                if result.eigenvalues is not None and len(result.eigenvalues) > n_components:
                                    result.eigenvalues = result.eigenvalues[:n_components]
                                if result.modes is not None and result.modes.shape[1] > n_components:
                                    result.modes = result.modes[:, :n_components]
                                if result.eigenvectors is not None and result.eigenvectors.shape[1] > n_components:
                                    result.eigenvectors = result.eigenvectors[:, :n_components]
                    except (ValueError, TypeError) as e:
                        # Skip problematic eigenvector data for now
                        pass

    # If we found cumvar, return the result object
    if result.cumvar is not None:
        return result

    # Fallback - return the original dataset list
    if dtype == 'dataset':
        return c_dslist
    else:
        return get_data_from_dtype(c_dslist, dtype)