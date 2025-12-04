"""
Visualization and mapping functions
"""
from .base import *

__all__ = ['volmap', 'rdf', 'pairdist', 'density', 'gist']


@super_dispatch()
def volmap(traj,
           mask,
           grid_spacing,
           size=None,
           center=None,
           buffer=3.0,
           centermask='*',
           radscale=1.36,
           peakcut=0.05,
           top=None,
           dtype='ndarray',
           frame_indices=None,
           options=""):
    '''(combined with cpptraj doc) Grid data as a volumetric map, similar to the
    volmap command in VMD. The density is calculated by treating each atom as a
    3-dimensional Gaussian function whose standard deviation is equal to the van der Waals radius

    Parameters
    ----------
    mask : {str, array-like}, default all atoms
        the atom selection from which to calculate the number density
    grid_spacing : tuple, grid spacing in X-, Y-, Z-dimensions, require
    size : {None, tuple}, default None
        if tuple, size must have length of 3
    center : {None, tuple}, default None
        if not None, center is tuple of (x, y, z) of center point
    buffer : float, default 3.0 Angstrom
        buffer distance (Angstrom), by which the edges of the grid should clear every atom
        of the centermask (or default mask if centermask is omitted) in every direction.
        The buffer is ignored if the center and size are specified.
    centermask : str
    radscale : float, default 1.36 (to match to VMD calculation)
        factor by which to scale radii (by devision)
    peakcut : float
    dtype : str, default 'ndarray'
        Note: To get all the output from cpptraj, it would be better to specify
        dtype='dict'

    Examples
    --------
    >>> import pytraj as pt
    >>> # load all frames to memory
    >>> traj = pt.datafiles.load_tz2_ortho()[:]
    >>> # do fitting and centering before perform volmap
    >>> traj = traj.superpose(mask=':1-13').center(':1-13 mass origin')
    >>> data = pt.volmap(traj, mask=':WAT@O', grid_spacing=(0.5, 0.5, 0.5), buffer=2.0, centermask='!:1-13', radscale=1.36)
    '''
    dummy_filename = 'dummy_fn.dat'

    if not isinstance(grid_spacing, tuple) or len(grid_spacing) != 3:
        raise ValueError('grid_spacing must be a tuple with length=3')

    if traj is None:
        raise ValueError('trajectory is required')

    grid_spacing_str = ' '.join([str(x) for x in grid_spacing])
    radscale_str = 'radscale ' + str(radscale)
    buffer_str = 'buffer ' + str(buffer)
    peakcut_str = 'peakcut ' + str(peakcut)
    centermask_str = 'centermask ' + centermask

    if isinstance(size, tuple):
        if len(size) != 3:
            raise ValueError('length of size must be 3')
    elif size is not None:
        raise ValueError(
            'size must be None or a tuple. Please check method doc')

    size_str = '' if size is None else 'size ' + ','.join(
        [str(x) for x in size])

    if size_str:
        # ignore buffer
        buffer_str = ''
        # add center
        if center is not None:
            center_str = 'center ' + ','.join([str(x) for x in center])
        else:
            center_str = ''
    else:
        center_str = ''

    command = ' '.join(
        (dummy_filename, grid_spacing_str, center_str, size_str, mask,
         radscale_str, buffer_str, centermask_str, peakcut_str, options))
    action_datasets, _ = do_action(traj, command, c_action.Action_Volmap)
    index = None
    for i, volume_ds in enumerate(action_datasets):
        if volume_ds.key.endswith("[totalvol]"):
            index = i
    if index is not None:
        action_datasets.remove_at(index)
    return get_data_from_dtype(action_datasets, dtype)


def rdf(traj=None,
        solvent_mask=':WAT@O',
        solute_mask='',
        maximum=10.,
        bin_spacing=0.5,
        image=True,
        density=0.033456,
        volume=False,
        center_solvent=False,
        center_solute=False,
        intramol=True,
        frame_indices=None,
        top=None,
        raw_rdf=False,
        dtype='tuple'):
    '''compute radial distribtion function. Doc was adapted lightly from cpptraj doc

    Returns
    -------
    a tuple of bin_centers, rdf values

    Parameters
    ----------
    traj : Trajectory-like, require
    solvent_mask : solvent mask, default None, required
    bin_spacing : float, default 0.5, optional
        bin spacing
    maximum : float, default 10., optional
        max bin value
    solute_mask: str, default None, optional
        if specified, calculate RDF of all atoms in solvent_mask to each atom in solute_mask
    image : bool, default True, optional
        if False, do not image distance
    density : float, default 0.033456 molecules / A^3, optional
    volume : determine density for normalization from average volume of input frames
    center_solvent : bool, default False, optional
        if True, calculate RDF from geometric center of atoms in solvent_mask to all atoms in solute_mask
    center_solute : bool, default False, optional
        if True, calculate RDF from geometric center of atoms in solute_mask to all atoms in solvent_mask
    intramol : bool, default True, optional
        if False, ignore intra-molecular distances
    frame_indices : array-like, default None, optional
    raw_rdf : bool, default False, optional
        if True, return the raw (non-normalized) RDF values
    dtype : str, default 'tuple'
        Return data type. 'tuple' returns (bin_centers, values), 'dict' returns
        {'bin_centers': ..., 'rdf': ...}, other types passed to get_data_from_dtype

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5, maximum=10.0, solute_mask=':WAT@O')
    >>> data[0]
    array([ 0.25,  0.75,  1.25, ...,  8.75,  9.25,  9.75])
    >>> data[1]
    array([ 0.        ,  0.        ,  0.        , ...,  0.95620052,
            0.95267934,  0.95135242])

    >>> # use array-like mask
    >>> atom_indices = pt.select(':WAT@O', traj.top)
    >>> data = pt.rdf(traj, solvent_mask=':WAT@O', bin_spacing=0.5, maximum=10.0, solute_mask=atom_indices)

    Notes
    -----
    - install ``pytraj`` and ``libcpptraj`` with openmp to speed up calculation
    - do not use this method with pytraj.pmap
    '''

    traj = get_fiterator(traj, frame_indices)

    # Convert masks to string format if needed
    if not isinstance(solvent_mask, str):
        solvent_mask = array_to_cpptraj_atommask(solvent_mask)

    if not isinstance(solute_mask, str) and solute_mask is not None:
        solute_mask = array_to_cpptraj_atommask(solute_mask)

    # Build command using CommandBuilder
    command = (CommandBuilder().add("pytraj_tmp_output.agr").add(
        str(bin_spacing)).add(str(maximum)).add(solvent_mask).add(
            solute_mask, condition=solute_mask
            is not None).add("noimage", condition=not image).add(
                "density", str(density), condition=density
                is not None).add("volume", condition=volume).add(
                    "center1", condition=center_solvent).add(
                        "center2", condition=center_solute).add(
                            "nointramol", condition=not intramol).add(
                                "rawrdf pytraj_tmp_output_raw.agr",
                                condition=raw_rdf).build())

    c_dslist, _ = do_action(traj, command, c_action.Action_Radial)
    # make a copy sine c_dslist[-1].values return view of its data
    # c_dslist will be freed
    values = np.array(c_dslist[-1].values)
    bin_centers = np.arange(bin_spacing / 2., maximum, bin_spacing)

    if dtype == 'tuple':
        # default behavior - return (bin_centers, values) tuple for backward compatibility
        return (bin_centers, values)
    elif dtype == 'dict':
        return {'bin_centers': bin_centers, 'rdf': values}
    else:
        # for 'ndarray', 'dataset', etc.
        return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def pairdist(traj,
             mask="*",
             mask2='',
             delta=0.1,
             dtype='ndarray',
             top=None,
             frame_indices=None):
    '''compute pair distribution function

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default all atoms
    mask2 : str, default ''
        second mask for pair distribution
    delta : float, default 0.1
        bin spacing
    dtype : str, default 'ndarray'
        dtype of return data
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.pairdist(traj)
    '''
    if traj is None:
        raise ValueError('trajectory is required')

    # Convert delta to float if it's a string and validate
    try:
        delta_val = float(delta)
        if delta_val <= 0:
            raise ValueError('delta must be positive')
    except (ValueError, TypeError):
        raise ValueError('delta must be a positive number')

    command = (CommandBuilder().add("mask",
                                    mask).add("mask2",
                                              mask2,
                                              condition=bool(mask2)).add(
                                                  "delta", str(delta)).build())

    action_datasets, _ = do_action(traj, command, c_action.Action_PairDist)
    return get_data_from_dtype(action_datasets, dtype=dtype)


def density(traj,
            mask='*',
            density_type='number',
            delta=0.25,
            direction='z',
            dtype='dict'):
    """Compute density (number, mass, charge, electron) along a coordinate

    Notes
    -----
    Syntax might be changed

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or list of str, default '*'
        required mask
    density_type : str, {'number', 'mass', 'charge', 'electron'}, default 'number'
    delta : float, default 0.25
        resolution (Angstrom)
    direction : str, default 'z'
    dtype : str, default 'dict'
        return data type. Please always using default value, others are for debugging.

    Returns
    -------
    out : dict of average density and std for each frame

    Examples
    --------

    >>> def func():
    ...     import pytraj as pt
    ...     fn = "data/DOPC.rst7"
    ...     tn = "data/DOPC.parm7"
    ...     traj = pt.load("data/DOPC.rst7", "data/DOPC.parm7")

    ...     delta = '0.25'
    ...     density_type = 'charge'
    ...     masks = [":PC@P31", ":PC@N31", ":PC@C2", ":PC | :OL | :OL2"]
    ...     density_dict = pt.density(traj, mask=masks, density_type=density_type, delta=delta)
    ...     return density_dict
    >>> density_dict = func() # doctest: +SKIP
    """

    assert density_type.lower() in {'number', 'mass', 'charge', 'electron'}, \
        f'{density_type} must be one of number, mass, charge, electron'

    if isinstance(mask, str):
        formatted_mask = f'"{mask}"'
    elif isinstance(mask, (list, tuple)):
        formatted_mask = ' '.join(f'"{m}"' for m in mask)
    else:
        raise ValueError("mask must be either string or list/tuple of string")

    command = f'delta {delta} {direction} {density_type} {formatted_mask}'
    action_datasets, _ = do_action(traj, command, c_action.Action_Density)

    result = get_data_from_dtype(action_datasets, dtype=dtype)
    if isinstance(result, dict):
        result.update({direction: action_datasets[0]._coord(dim=0)})

    return result


@super_dispatch()
def gist(traj,
         solute_mask='',
         solvent_mask='',
         gridcenter=(0, 0, 0),
         griddim=(10, 10, 10),
         gridspacn=0.5,
         dtype='dataset',
         top=None,
         frame_indices=None):
    """compute GIST analysis

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str, optional
        solute mask
    solvent_mask : str, optional
        solvent mask
    gridcenter : tuple of floats, default (0, 0, 0)
        grid center coordinates
    griddim : tuple of ints, default (10, 10, 10)
        grid dimensions
    gridspacn : float, default 0.5
        grid spacing
    dtype : str, default 'dataset'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : DatasetList
    """
    cx, cy, cz = gridcenter
    nx, ny, nz = griddim

    command = f"gridcenter {cx} {cy} {cz} griddim {nx} {ny} {nz} gridspacn {gridspacn}"
    if solute_mask:
        command = f"solute {solute_mask} " + command
    if solvent_mask:
        command = f"solvent {solvent_mask} " + command

    action_datasets, _ = do_action(traj, command, c_action.Action_Density)
    return get_data_from_dtype(action_datasets, dtype=dtype)


def _grid(traj,
          size=(1., 1., 1.),
          center=(0, 0, 0),
          mask='*',
          buffer=2.0,
          max_peak=0,
          npoints=40000,
          invert=False):
    """compute grid data (density) around center or center of mask

    Parameters
    ----------
    traj : Trajectory-like
    size : tuple of floats, default (1., 1., 1.)
        grid size
    center : tuple of floats, default (0, 0, 0)
        center coordinates
    mask : str, default '*'
        if given and center=(0, 0, 0), use mask center as center
    buffer : float, default 2.0
        buffer distance
    max_peak : int, default 0
        number of high density peaks to find
    npoints : int, default 40000
        number of grid points
    invert : bool, default False
        if True, invert density

    Returns
    -------
    out : ndarray, coordinates of density peak
    """
    x_size, y_size, z_size = size
    cx, cy, cz = center

    if cx == 0. and cy == 0. and cz == 0.:
        # use mask center
        command = f"grid {x_size} {y_size} {z_size} {mask} max {max_peak} buffer {buffer} npoints {npoints}"
    else:
        command = f"grid {x_size} {y_size} {z_size} {cx} {cy} {cz} max {max_peak} buffer {buffer} npoints {npoints}"

    if invert:
        command += " invert"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_Grid()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return c_dslist[-1].values.reshape(-1, 3)
