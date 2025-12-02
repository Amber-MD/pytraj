"""
Visualization and mapping functions
"""
from .base import *

__all__ = [
    'volmap', 'rdf', 'pairdist', 'density', 'gist'
]


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

    assert isinstance(grid_spacing, tuple) and len(grid_spacing) == 3, 'grid_spacing must be a tuple with length=3'

    grid_spacing_str = ' '.join([str(x) for x in grid_spacing])
    radscale_str = 'radscale ' + str(radscale)
    buffer_str = 'buffer ' + str(buffer)
    peakcut_str = 'peakcut ' + str(peakcut)
    centermask_str = 'centermask ' + centermask

    if isinstance(size, tuple):
        assert len(size) == 3, 'length of size must be 3'
    elif size is not None:
        raise ValueError('size must be None or a tuple. Please check method doc')

    size_str = '' if size is None else 'size ' + ','.join([str(x) for x in size])

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

    command = ' '.join((dummy_filename, grid_spacing_str, center_str, size_str, mask,
                        radscale_str, buffer_str, centermask_str, peakcut_str, options))
    action_datasets, _ = do_action(traj, command, c_action.Action_Volmap)
    index = None
    for i, volume_ds in enumerate(action_datasets):
        if volume_ds.key.endswith("[totalvol]"):
            index = i
    if index is not None:
        action_datasets._pop(index)
    return get_data_from_dtype(action_datasets, dtype)


def rdf(traj=None,
        solute_mask='',
        solvent_mask='',
        bin_spacing=0.5,
        bin_max=10.0,
        volume=False,
        density=0.033456,
        center1=False,
        center2=False,
        intramol=True,
        image=True,
        dtype='ndarray',
        top=None,
        frame_indices=None):
    """compute radial distribution functions

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask : str
        solute mask
    solvent_mask : str
        solvent mask
    bin_spacing : float, default 0.5
        histogram bin spacing
    bin_max : float, default 10.0
        histogram max distance
    volume : bool, default False
        if True, normalize by volume
    density : float, default 0.033456
        number density (molecules per cubic Angstrom)
    center1 : bool, default False
        if True, only calculate to center of first mask
    center2 : bool, default False
        if True, only calculate to center of second mask
    intramol : bool, default True
        if True, include intramolecular interactions
    image : bool, default True
        if True, use nearest image
    dtype : str, default 'ndarray'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : ndarray, shape (n_bins, 2)
    """
    command = f"{solute_mask} {solvent_mask} {bin_spacing} {bin_max}"

    if volume:
        command += " volume"
    command += f" density {density}"
    if center1:
        command += " center1"
    if center2:
        command += " center2"
    if not intramol:
        command += " nointramol"
    if not image:
        command += " noimage"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_Radial()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def pairdist(traj,
             mask='',
             ref_mask='',
             dtype='dataset',
             top=None,
             frame_indices=None):
    """compute pairwise distances

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom selection
    ref_mask : str, optional
        reference mask to compute distances to
    dtype : str, default 'dataset'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : DatasetList or ndarray
    """
    command = mask
    if ref_mask:
        command += f" ref {ref_mask}"

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_Pairwise()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def density(traj,
            mask='*',
            density_type='number',
            delta=0.25,
            direction='z',
            dtype='dict',
            top=None):
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

    c_dslist = CpptrajDatasetList()
    action = c_action.Action_Gist()
    action.read_input(command, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


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