from __future__ import absolute_import
import numpy as np

from pytraj.trajectory import Trajectory
from pytraj.trajectory_iterator import TrajectoryIterator
from .get_common_objects import (get_topology, get_data_from_dtype, get_list_of_commands,
                                 get_matrix_from_dataset, get_reference, get_fiterator,
                                 super_dispatch, get_fi_with_dslist)
from .utils import ensure_not_none_or_string
from .utils import is_int
from .utils.context import tempfolder
from .utils.convert import array_to_cpptraj_atommask
from .externals.six import string_types
from .datasets.c_datasetlist import DatasetList as CpptrajDatasetList
from .datasetlist import DatasetList
from .shared_methods import iterframe_master
from .decorators import register_pmap, register_openmp
from .c_action import c_action
from .c_analysis import c_analysis
from .c_action.actionlist import ActionList
from .utils.convert import array2d_to_cpptraj_maskgroup
from .topology import Topology

list_of_calc = ['calc_distance',
                'calc_dihedral',
                'calc_radgyr',
                'calc_angle',
                'calc_surf',
                'calc_molsurf',
                'calc_volume',
                'calc_matrix',
                'calc_jcoupling',
                'calc_watershell',
                'calc_vector',
                'calc_multivector',
                'calc_volmap',
                'calc_rdf',
                'calc_pairdist',
                'calc_multidihedral',
                'calc_atomicfluct',
                'calc_center_of_mass',
                'calc_center_of_geometry',
                'calc_pairwise_rmsd',
                'calc_grid',
                'calc_atomiccorr',
                'calc_bfactors',
                'calc_diffusion',
                'calc_distance_rmsd',
                'calc_mindist',
                'calc_pairwise_distance',
                'calc_rmsd_nofit',
                'calc_rotation_matrix',
                'calc_pca', ]

list_of_calc_short = [word.replace('calc_', '')
                      for word in list_of_calc if not word in ['calc_matrix', ]]

list_of_do = ['translate', 'rotate', 'autoimage', 'scale']

list_of_get = ['get_average_frame', 'get_velocity']

list_of_the_rest = ['rmsd', 'align_principal_axis', 'principal_axes', 'closest',
                    'transform', 'native_contacts', 'set_dihedral',
                    'check_structure', 'mean_structure', 'lowestcurve',
                    'make_structure', 'replicate_cell', 'pucker', 'rmsd_perres',
                    'randomize_ions', 'velocityautocorr',
                    'timecorr', 'search_neighbors',
                    'xcorr', 'acorr',
                    'projection',
                    'superpose', 'strip',
                    'center', 'wavelet'
                    ]

__all__ = list(set(list_of_do + list_of_calc_short + list_of_get + list_of_the_rest))


def _2darray_to_atommask_groups(seq):
    '''

    >>> list(_2darray_to_atommask_groups([[0, 3], [4, 7]]))
    ['@1 @4', '@5 @8']
    '''
    for arr in seq:
        # example: arr = [0, 3]; turns ot '@1 @4'
        yield '@' + str(arr[0] + 1) + ' @' + str(arr[1] + 1)


def _assert_mutable(trajiter):
    """make sure the input is not TrajectoryIterator
    """
    if isinstance(trajiter, TrajectoryIterator):
        raise ValueError(
            "This analysis does not support immutable object. Use `pytraj.Trajectory`")


@register_pmap
def distance(traj=None,
             mask="",
             frame_indices=None,
             dtype='ndarray',
             top=None,
             image=False,
             n_frames=None):
    # TODO: add image, noe, ...
    """compute distance between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or a list of string or a 2D array-like of integers
    frame_indices : array-like, optional, default None
    dtype : return type, default 'ndarray'
    top : Topology, optional
    image : bool, default False
    n_frames : int, optional, default None
        only need to provide n_frames if ``traj`` does not have this info

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Notes
    -----
    Be careful with Topology. If your topology has Box info but your traj does not, you
    would get weird output ([0.0, ...]). Make sure to use `image=False` in this method or
    set_nobox for Topology.


    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate distance for two atoms, using amber mask
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> dist = pt.distance(traj, '@1 @3')

    >>> # calculate distance for two groups of atoms, using amber mask
    >>> dist = pt.distance(traj, '@1,37,8 @2,4,6')

    >>> # calculate distance between two residues, using amber mask
    >>> dist = pt.distance(traj, ':1 :10')

    >>> # calculate multiple distances between two residues, using amber mask
    >>> # distance between residue 1 and 10, distance between residue 3 and 20
    >>> # (when using atom string mask, index starts from 1)
    >>> dist = pt.distance(traj, [':1 :10', ':3 :20'])

    >>> # calculate distance for a series of atoms, using array for atom mask
    >>> # distance between atom 1 and 5, distance between atom 4 and 10 (index starts from 0)
    >>> dist = pt.distance(traj, [[1, 5], [4, 10]])
    """
    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    top_ = get_topology(traj, top)
    _noimage = 'noimage' if not image else ''

    cm_arr = np.asarray(command)

    if 'int' in cm_arr.dtype.name:

        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(cm_arr)

        n_frames = traj.n_frames if n_frames is None else n_frames

        arr = np.empty([n_frames, len(int_2darr)])

        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._distance(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            dslist = DatasetList({'distance': arr})
            return get_data_from_dtype(dslist, dtype)

    elif isinstance(command, (list, tuple, string_types, np.ndarray)):
        # create a list
        if not isinstance(command, np.ndarray):
            list_of_commands = get_list_of_commands(command)
        else:
            list_of_commands = command

        c_dslist = CpptrajDatasetList()
        actlist = ActionList()

        for cm in list_of_commands:
            if not image:
                cm = ' '.join((cm, _noimage))
            actlist.add(c_action.Action_Distance(),
                        cm,
                        top_,
                        dslist=c_dslist)

        actlist.compute(traj)
        return get_data_from_dtype(c_dslist, dtype)

    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")

calc_distance = distance


def pairwise_distance(traj=None,
                      mask_1='',
                      mask_2='',
                      top=None,
                      dtype='ndarray',
                      frame_indices=None):
    '''compute pairwise distance between atoms in mask_1 and atoms in mask_2

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask_1: string or 1D array-like
    mask_2: string or 1D array-like
    ...

    Returns
    -------
    out_1 : numpy array, shape (n_frames, n_atom_1, n_atom_2)
    out_2 : atom pairs, shape=(n_atom_1, n_atom_2, 2)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> mat = pt.calc_pairwise_distance(traj, '@CA', '@CB')

    Notes
    -----
    This method is only fast for small number of atoms.
    '''
    from itertools import product

    top_ = get_topology(traj, top)
    indices_1 = top_.select(mask_1) if isinstance(mask_1,
                                                  string_types) else mask_1
    indices_2 = top_.select(mask_2) if isinstance(mask_2,
                                                  string_types) else mask_2
    arr = np.array(list(product(indices_1, indices_2)))
    mat = calc_distance(traj,
                        mask=arr,
                        dtype=dtype,
                        top=top_,
                        frame_indices=frame_indices)
    mat = mat.T
    return (mat.reshape(mat.shape[0], len(indices_1), len(indices_2)),
            arr.reshape(
                len(indices_1), len(indices_2), 2))

calc_pairwise_distance = pairwise_distance


@register_pmap
def angle(traj=None,
          mask="",
          top=None,
          dtype='ndarray',
          frame_indices=None,
          *args,
          **kwargs):
    """compute angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # calculate angle for three atoms, using amber mask
    >>> pt.angle(traj, '@1 @3 @10')
    array([  98.06193365,   95.75979717,  105.59626997,  107.64079091,
             94.93516228,  102.06028369,  109.3479057 ,  110.11532478,
            101.86104127,  110.48992512])

    >>> # calculate angle for three groups of atoms, using amber mask
    >>> angles = pt.angle(traj, '@1,37,8 @2,4,6 @5,20')

    >>> # calculate angle between three residues, using amber mask
    >>> angles = pt.angle(traj, ':1 :10 :20')

    >>> # calculate multiple angles between three residues, using amber mask
    >>> # angle between residue 1, 10, 20, angle between residue 3, 20, 30
    >>> # (when using atom string mask, index starts from 1)
    >>> angles = pt.angle(traj, [':1 :10 :20', ':3 :20 :30'])

    >>> # calculate angle for a series of atoms, using array for atom mask
    >>> # angle between atom 1, 5, 8, distance between atom 4, 10, 20 (index starts from 0)
    >>> angles = pt.angle(traj, [[1, 5, 8], [4, 10, 20]])
    """
    c_dslist = CpptrajDatasetList()
    act = c_action.Action_Angle()

    command = mask

    ensure_not_none_or_string(traj)

    traj = get_fiterator(traj, frame_indices)
    top_ = get_topology(traj, top)
    cm_arr = np.asarray(command)

    if 'int' not in cm_arr.dtype.name:
        if isinstance(command, string_types):
            # need to remove 'n_frames' keyword since Action._master does not use
            # it
            try:
                del kwargs['n_frames']
            except KeyError:
                pass
            # cpptraj mask for action
            act(command, traj, top=top_, dslist=c_dslist, *args, **kwargs)
            return get_data_from_dtype(c_dslist, dtype)
        elif isinstance(command, (list, tuple, np.ndarray)):
            list_of_commands = command
            c_dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add(c_action.Action_Angle(),
                            cm,
                            top_,
                            dslist=c_dslist,
                            *args,
                            **kwargs)
            actlist.compute(traj)
            return get_data_from_dtype(c_dslist, dtype)

    else:
        # ndarray of integer
        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(int_2darr)

        if int_2darr.shape[1] != 3:
            raise ValueError("require int-array with shape=(n_atoms, 3)")

        n_frames = kwargs.get('n_frames')
        n_frames = traj.n_frames if n_frames is None else n_frames

        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._angle(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = DatasetList({'angle': arr})
            return get_data_from_dtype(py_dslist, dtype)

calc_angle = angle


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    '''compute dihedral within a single residue. For internal use only.

    Parameters
    ----------
    traj : Trajectory-like
    mask : tuple of strings
    resid : str, resid
    dtype

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.all_actions import _dihedral_res
    >>> traj = pt.datafiles.load_tz2()
    >>> data = _dihedral_res(traj, mask=('N', 'CA', 'C', 'O'), resid=0)
    >>> # use string for resid
    >>> data = _dihedral_res(traj, mask=('N', 'CA', 'C', 'O'), resid='1')
    '''

    if is_int(resid):
        resid = str(resid + 1)
    else:
        resid = resid
    m = ' :%s@' % resid
    command = m + m.join(mask)
    return calc_dihedral(traj=traj, mask=command, top=top, dtype=dtype)


@register_pmap
def dihedral(traj=None,
             mask="",
             top=None,
             dtype='ndarray',
             frame_indices=None,
             *args,
             **kwargs):
    """compute dihedral angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Returns
    -------
    1D ndarray if mask is a string
    2D ndarray, shape (n_atom_pairs, n_frames) if mask is a list of strings or an array

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # calculate dihedral angle for four atoms, using amber mask
    >>> pt.dihedral(traj, '@1 @3 @10 @20')
    array([ 23.32244362,  23.5386922 ,  14.26831569,  14.58865946,
            23.98675475,  26.18419185,   6.06982926,  13.57158505,
            16.59013076,  29.99131573])

    >>> # calculate dihedral angle for four groups of atoms, using amber mask
    >>> data = pt.dihedral(traj, '@1,37,8 @2,4,6 @5,20 @21,22')

    >>> # calculate dihedral angle for four residues, using amber mask
    >>> data = pt.dihedral(traj, ':1 :10 :20 :22')

    >>> # calculate multiple dihedral angles for four residues, using amber mask
    >>> # angle for residue 1, 10, 20, 30; angle between residue 3, 20, 30, 40
    >>> # (when using atom string mask, index starts from 1)
    >>> pt.dihedral(traj, [':1 :10 :20 :30', ':3 :20 :30 :40'])
    array([[-166.81829116, -160.19029669, -158.56560062, ..., -169.10064772,
            -158.81655586, -165.28873555],
           [  -0.60994639,    0.78261235,    1.86394369, ...,    1.21170489,
              -1.16762607,   -3.08412049]])
    >>> # calculate dihedral angle for a series of atoms, using array for atom mask
    >>> # dihedral angle for atom 1, 5, 8, 10, dihedral for atom 4, 10, 20, 30 (index starts from 0)
    >>> data = pt.dihedral(traj, [[1, 5, 8, 10], [4, 10, 20, 30]])
    """
    act = c_action.Action_Dihedral()
    c_dslist = CpptrajDatasetList()

    ensure_not_none_or_string(traj)
    command = mask

    traj = get_fiterator(traj, frame_indices)
    top_ = get_topology(traj, top)
    cm_arr = np.asarray(command)

    if 'int' not in cm_arr.dtype.name:
        if isinstance(command, string_types):
            # need to remove 'n_frames' keyword since Action._master does not use
            # it
            try:
                del kwargs['n_frames']
            except:
                pass
            # cpptraj mask for action
            act(command, traj, top=top_, dslist=c_dslist)
            return get_data_from_dtype(c_dslist, dtype)

        else:
            # assume isinstance(command, (list, tuple, np.ndarray)):
            list_of_commands = command
            c_dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add(c_action.Action_Dihedral(),
                            cm,
                            top_,
                            dslist=c_dslist,
                            *args,
                            **kwargs)
            actlist.compute(traj)
            return get_data_from_dtype(c_dslist, dtype)
    else:
        # ndarray
        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(int_2darr)

        if int_2darr.shape[1] != 4:
            raise ValueError("require int-array with shape=(n_atoms, 4)")

        n_frames = kwargs.get('n_frames')
        n_frames = traj.n_frames if n_frames is None else n_frames
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._dihedral(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = DatasetList({'dihedral': arr})
            return get_data_from_dtype(py_dslist, dtype)

calc_dihedral = dihedral


@register_pmap
def mindist(traj=None,
            command="",
            top=None,
            dtype='ndarray',
            frame_indices=None):
    '''compute mindist

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.mindist(traj, '@CA @H')
    '''

    traj = get_fiterator(traj, frame_indices)
    act = c_action.Action_NativeContacts()
    c_dslist = CpptrajDatasetList()

    if not isinstance(command, string_types):
        command = array2d_to_cpptraj_maskgroup(command)
    command_ = "mindist " + command
    traj = get_fiterator(traj, frame_indices)
    top_ = get_topology(traj, top)
    act(command_, traj, top=top_, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype=dtype)[-1]

calc_mindist = mindist


@super_dispatch()
def diffusion(traj,
              mask="",
              tstep=1.0,
              individual=False,
              top=None,
              dtype='dataset',
              frame_indices=None):
    '''compute diffusion

    Parameters
    ----------
    traj : Trajectory-like or iterator
    mask : str, default '' (all atoms)
    tstep : float, time step between frames, default 1.0 ps
    individual : bool, default False
    top : Topology, optional
    dtype : str, default 'dataset'
    frame_indices : array or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.diffusion(traj, dtype='dict')
    >>> data['X']
    array([ 0.        ,  0.87027302,  1.64626022,  2.26262651,  2.98068114,
            3.57075535,  4.07030655,  4.71894512,  5.42302306,  6.01310377])
    '''
    act = c_action.Action_Diffusion()
    c_dslist = CpptrajDatasetList()

    _tsep = 'time ' + str(tstep)
    _individual = 'individual' if individual else ''

    # add 'df' as label
    label = 'df'
    command = ' '.join((mask, label, _tsep, _individual))

    # normally we just need
    # act(command, traj, top=top_, dslist=c_dslist)
    # but cpptraj need correct frame idx

    act.read_input(command, top=top, dslist=c_dslist)
    act.setup(top)
    for frame in traj:
        # do not need mass
        act.compute(frame)
    act.post_process()

    # make the label nicer
    for d in c_dslist:
        d.key = d.key.replace('[', '').replace(']', '').replace(label, '')

    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_diffusion = diffusion


@register_pmap
@register_openmp
def watershell(traj=None,
               solute_mask='',
               solvent_mask=':WAT',
               lower=3.4,
               upper=5.0,
               image=True,
               dtype='dataset',
               frame_indices=None,
               top=None):
    """(adapted from cpptraj doc): Calculate numbers of waters in 1st and 2nd solvation shells
    (defined by <lower cut> (default 3.4 Ang.) and <upper cut> (default 5.0 Ang.)

    Parameters
    ----------
    traj : Trajectory-like
    solute_mask: solute mask
    solvent_mask: solvent mask
    lower : double, default 3.4
        lower cut distance
    upper : double, default 5.0
        upper cut distance
    image : bool, defaul True
        do autoimage if True
    dtype : return type, defaul 'dataset'
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.watershell(traj, solute_mask='!:WAT')
    >>> data = pt.watershell(traj, solute_mask='!:WAT', lower=5.0, upper=10.)
    """
    _solutemask = solute_mask
    top_ = get_topology(traj, top)
    fi = get_fiterator(traj, frame_indices)

    c_dslist = CpptrajDatasetList()

    act = c_action.Action_Watershell()

    if _solutemask in [None, '']:
        raise ValueError('must provide solute mask')

    _solventmask = solvent_mask if solvent_mask is not None else ''
    _noimage = 'noimage' if not image else ''

    _lower = 'lower ' + str(lower)
    _upper = 'upper ' + str(upper)
    command = ' '.join((_solutemask, _lower, _upper, _noimage, _solventmask))

    act(command, fi, top=top_, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_watershell = watershell


@super_dispatch()
def calc_matrix(traj=None,
                mask="",
                top=None,
                dtype='ndarray',
                frame_indices=None):
    '''compute different type of matrices

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, type of matrix and atom mask
    top : Topology, optional
    dtype: return data type
    frame_indices : {None, array-like}
        if not None, perform calculation for given frame indices

    Notes
    -----
    If user wants to use specify matrix's method name, see also ``pytraj.matrix``

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_trpcage()
    >>> mat = pt.calc_matrix(traj, 'covar @CA')
    >>> # this is equal to
    >>> mat2 = pt.matrix.covar(traj, '@CA')
    >>> import numpy as np
    >>> np.testing.assert_equal(mat, mat2)
    '''
    act = c_action.Action_Matrix()

    c_dslist = CpptrajDatasetList()
    act(mask, traj, top=top, dslist=c_dslist)
    act.post_process()
    return get_data_from_dtype(c_dslist, dtype)


@register_pmap
@super_dispatch()
def radgyr(traj=None,
           mask="",
           top=None,
           nomax=True,
           frame_indices=None,
           dtype='ndarray'):
    '''compute radius of gyration

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr(traj, '@CA')
    >>> data = pt.radgyr(traj, '!:WAT', nomax=False)
    >>> data = pt.radgyr(traj, '@CA', frame_indices=[2, 4, 6])
    '''
    nomax_ = 'nomax' if nomax else ""
    command = " ".join((mask, nomax_))

    act = c_action.Action_Radgyr()

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype)

calc_radgyr = radgyr


@register_pmap
@super_dispatch()
def surf(traj=None,
         mask="",
         dtype='ndarray',
         frame_indices=None,
         top=None):
    '''calc surf (LCPO method)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.surf(traj, '@CA')
    '''
    act = c_action.Action_Surf()

    command = mask

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype)

calc_surf = surf


@register_pmap
@super_dispatch()
def molsurf(traj=None,
            mask="",
            probe=1.4,
            offset=0.0,
            dtype='ndarray',
            frame_indices=None,
            top=None):
    '''calc molsurf

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> pt.molsurf(traj, '@CA')
    array([ 458.51409637,  459.64784573,  456.54690793,  467.72939574,
            462.45908781,  458.70327554,  454.40514806,  455.15015576,
            468.70566447,  456.0058624 ])
    >>> pt.molsurf(traj, '!:WAT')
    array([ 1079.1395679 ,  1090.79759341,  1069.65127413,  1096.0810919 ,
            1091.65862234,  1091.68906298,  1085.53105392,  1069.22510187,
            1079.70803583,  1075.8151414 ])
    '''
    _probe = 'probe ' + str(probe)
    offset_ = 'offset ' + str(offset) if offset != 0. else ''
    command = ' '.join((mask, _probe, offset_))

    act = c_action.Action_Molsurf()

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype)

calc_molsurf = molsurf


@register_pmap
@super_dispatch()
def rotation_matrix(traj=None,
                    mask="",
                    ref=0,
                    mass=False,
                    frame_indices=None,
                    top=None,
                    with_rmsd=False):
    '''compute rotation matrix with/without rmsd

    Parameters
    ----------
    traj : Trajectory-like
    ref : {int, Frame}, default 0 (first Frame)
        reference
    mask : str, default all atoms
    mass : bool, default False
        if True, rmsfit with mass
    frame_indices : {None, array-like}
        if not None, compute for given indices
    top : Topology, optional
    with_rmsd : bool, default False
        - if False, return only rotation matrix.
        - if True, return rotation matrix and rmsd values

    Returns
    -------
    out : if with_rmsd=False, return numpy array, shape (n_frames, 3, 3)
          if with_rmsd=True, return a tuple (mat, rmsd)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> mat = pt.calc_rotation_matrix(traj, mask='@CA')
    >>> mat.shape
    (101, 3, 3)
    '''
    c_dslist = CpptrajDatasetList()
    mass_ = 'mass' if mass else ''

    command = ' '.join(('tmp', mask, 'savematrices', mass_))

    act = c_action.Action_Rmsd()
    act(command, [ref, traj], top=top, dslist=c_dslist)
    mat = c_dslist[-1].values
    # exclude data for reference
    if with_rmsd:
        return mat[1:], np.array(c_dslist[0].values[1:])
    else:
        return mat[1:]

calc_rotation_matrix = rotation_matrix


@super_dispatch()
def volume(traj=None,
           mask="",
           top=None,
           dtype='ndarray',
           frame_indices=None):
    '''compute volume

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> vol = pt.calc_volume(traj, '@CA')
    '''
    command = mask

    act = c_action.Action_Volume()

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype)

calc_volume = volume


@super_dispatch()
def multivector(traj,
                resrange,
                names,
                top=None,
                dtype='dataset',
                frame_indices=None):
    '''

    Parameters
    ----------
    traj : Trajectory-like
    resrange : str, residue range
    names : {str, tuple of str}
    top : Topology, optional
    dtype : str, default 'dataset'
    frame_indices : {None, 1D array-like}, optional, default None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> vecs = pt.multivector(traj, resrange='1-5', names=('C', 'N'))
    >>> vecs = pt.multivector(traj, resrange='1-5', names='C N')
    '''
    act = c_action.Action_MultiVector()

    _resrange = 'resrange ' + resrange
    if 'name1' in names or 'name2' in names:
        # cpptraj style
        _names = names
    else:
        if isinstance(names, string_types):
            name1, name2 = names.split()
        else:
            # try to unpack
            name1, name2 = names
        _names = ' '.join(('name1', name1, 'name2', name2))
    command = ' '.join((_resrange, _names))

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype)

calc_multivector = multivector


@register_pmap
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
           frame_indices=None):
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

    grid_spacing_ = ' '.join([str(x) for x in grid_spacing])
    radscale_ = 'radscale ' + str(radscale)
    buffer_ = 'buffer ' + str(buffer)
    peakcut_ = 'peakcut ' + str(peakcut)
    centermask_ = 'centermask ' + centermask

    if isinstance(size, tuple):
        assert len(size) == 3, 'lenghth of size must be 3'
    elif size is not None:
        raise ValueError('size must be None or a tuple. Please check method doc')

    size_ = '' if size is None else 'size ' + ','.join([str(x) for x in size])

    if size_:
        # ignore buffer
        buffer_ = ''
        # add center
        if center is not None:
            center_ = 'center ' + ','.join([str(x) for x in center])
        else:
            center_ = ''
    else:
        center_ = ''

    command = ' '.join((dummy_filename, grid_spacing_, center_, size_, mask, radscale_, buffer_,
                        centermask_, peakcut_))

    act = c_action.Action_Volmap()

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=c_dslist)
    act.post_process()
    return get_data_from_dtype(c_dslist, dtype)


calc_volmap = volmap


@register_openmp
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
        top=None):
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
    top_ = get_topology(traj, top)

    act = c_action.Action_Radial()

    if not isinstance(solvent_mask, string_types):
        solvent_mask = array_to_cpptraj_atommask(solvent_mask)

    if not isinstance(solute_mask, string_types) and solute_mask is not None:
        solute_mask = array_to_cpptraj_atommask(solute_mask)

    spacing_ = str(bin_spacing)
    _maximum = str(maximum)
    _solventmask = solvent_mask
    _solutemask = solute_mask
    _noimage = 'noimage' if not image else ''
    _density = 'density ' + str(density) if density is not None else ''
    volume_ = 'volume' if volume else ''
    _center1 = 'center1' if center_solvent else ''
    _center2 = 'center2' if center_solute else ''
    _nointramol = 'nointramol' if not intramol else ''

    # order does matters
    # the order between _solventmask and _solutemask is swapped compared
    # to cpptraj's doc (to get correct result)
    command = ' '.join(
        ("pytraj_tmp_output.agr", spacing_, _maximum, _solventmask, _solutemask,
         _noimage, _density, volume_, _center1, _center2, _nointramol))

    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top_, dslist=c_dslist)
    act.post_process()

    # make a copy sine c_dslist[-1].values return view of its data
    # c_dslist will be freed
    values = np.array(c_dslist[-1].values)
    # return (bin_centers, values)
    return (np.arange(bin_spacing / 2., maximum, bin_spacing), values)

calc_rdf = rdf


@super_dispatch()
def calc_pairdist(traj,
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
    delta : float, default 0.1
    dtype : str, default 'ndarray'
        dtype of return data
    top : Topology, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.calc_pairdist(traj)
    '''
    c_dslist = CpptrajDatasetList()
    act = c_action.Action_PairDist()

    mask_ = 'mask ' + mask
    mask2_ = 'mask2 ' + str(mask2) if mask2 else ''
    delta_ = 'delta ' + str(delta)
    command = ' '.join((mask_, mask2_, delta_))

    act(command, traj, top=top, dslist=c_dslist)
    act.post_process()

    return get_data_from_dtype(c_dslist, dtype=dtype)


pairdist = calc_pairdist


@super_dispatch()
def jcoupling(traj=None,
              mask="",
              top=None,
              kfile=None,
              dtype='dataset',
              frame_indices=None):
    """compute j-coupling

    Parameters
    ----------
    traj : any things that make `frame_iter_master` returning Frame object
    command : str, default ""
        cpptraj's command/mask
    kfile : str, default None, optional
        Dir for Karplus file. If "None", use $AMBERHOME dir
    dtype : str, {'dataset', ...}, default 'dataset'

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.calc_jcoupling(traj, ':1-12', kfile='data/Karplus.txt')
    """
    command = mask

    act = c_action.Action_Jcoupling()
    c_dslist = CpptrajDatasetList()

    if kfile is not None:
        command += " kfile %s" % kfile
    act(command, traj, dslist=c_dslist, top=top)
    return get_data_from_dtype(c_dslist, dtype)

calc_jcoupling = jcoupling


def translate(traj=None, command="", frame_indices=None, top=None):
    '''translate coordinate

    Examples
    --------
    >>> import pytraj as pt
    >>> # load to mutable trajectory by `load` method
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.load(fn, tn)
    >>> # do transform traj and return itself
    >>> traj = pt.translate(traj, '@CA x 120.')
    '''
    _assert_mutable(traj)

    top_ = get_topology(traj, top)
    fi = get_fiterator(traj, frame_indices)

    c_action.Action_Translate()(command, fi, top=top_)


do_translation = translate


def do_scaling(traj=None, command="", frame_indices=None, top=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> # load to mutable trajectory by `load` method
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.load(fn, tn)
    >>> # do transform traj and return itself
    >>> traj = pt.scale(traj, '@CA x 1.2')
    '''
    _assert_mutable(traj)
    top_ = get_topology(traj, top)
    fi = get_fiterator(traj, frame_indices)
    c_action.Action_Scale()(command, fi, top=top_)


scale = do_scaling


def rotate(traj=None, command="", frame_indices=None, top=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> # load to mutable trajectory by `load` method
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.load(fn, tn)
    >>> # do transform traj and return itself
    >>> traj = pt.rotate(traj, 'x 90')

    Notes
    -----
    ``rotate`` is an alias of ``do_rotation``
    '''
    top_ = get_topology(traj, top)
    _assert_mutable(traj)
    fi = get_fiterator(traj, frame_indices)
    c_action.Action_Rotate()(command, fi, top=top_)


do_rotation = rotate


@super_dispatch()
def do_autoimage(traj,
                 mask="",
                 frame_indices=None,
                 top=None):
    '''perform autoimage and return the coordinate-updated traj

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()[:]
    >>> traj = pt.autoimage(traj)
    '''
    _assert_mutable(traj)
    c_action.Action_AutoImage()(mask, traj, top=top)
    return traj


autoimage = do_autoimage


@register_pmap
def mean_structure(traj,
                   mask='',
                   frame_indices=None,
                   dtype='frame',
                   autoimage=False,
                   rmsfit=None,
                   top=None):
    '''get mean structure for a given mask and given frame_indices

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : None or str, default None (all atoms)
    frame_indices : iterable that produces integer, default None, optional
        frame indices
    dtype: str, {'frame', 'traj'}, default 'frame'
        return type, either Frame (does not have Topology information) or 'traj'
    autoimage : bool, default False
        if True, performa autoimage
    rmsfit : object, {Frame, int, tuple, None}, default None
        if rmsfit is not None, perform rms fit to reference.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # get average structure from whole traj, all atoms
    >>> frame = pt.mean_structure(traj)

    >>> # get average structure from several frames, '@CA' atoms"
    >>> frame = pt.mean_structure(traj, '@CA', frame_indices=range(2, 8, 2))

    >>> # get average structure but do autoimage and rmsfit to 1st Frame
    >>> frame = pt.mean_structure(traj(autoimage=True, rmsfit=0))

    >>> # get average structure but do autoimage and rmsfit to 1st Frame.
    >>> frame = pt.mean_structure(traj(autoimage=True, rmsfit=0, frame_indices=[0, 5, 6]))

    Notes
    -----
    if autoimage=True and having rmsfit, perform autoimage first and do rmsfit
    '''
    # note: we not yet use @super_dispatch due to extra 'rmsfit'
    # TODO: do it.
    top_ = get_topology(traj, top)
    try:
        fi = traj.iterframe(autoimage=autoimage,
                            rmsfit=rmsfit,
                            frame_indices=frame_indices)
    except AttributeError:
        fi = get_fiterator(traj, frame_indices)

    c_dslist = CpptrajDatasetList()
    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)
    else:
        mask = mask

    # add "crdset s1" to trick cpptraj dump coords to DatSetList
    command = mask + " crdset s1"

    act = c_action.Action_Average()
    act(command, fi, top_, dslist=c_dslist)

    # need to call this method so cpptraj will write
    act.post_process()

    frame = c_dslist[0].get_frame()
    if dtype.lower() == 'frame':
        return frame
    elif dtype.lower() in ['traj', 'trajectory']:
        new_top = top_ if mask is '' else top_[mask]
        return Trajectory(xyz=frame.xyz.reshape(1, frame.n_atoms, 3).copy(),
                          top=new_top)
    else:
        raise ValueError('dtype must be frame or trajectory')


get_average_frame = mean_structure


@register_pmap
def get_velocity(traj, mask=None, frame_indices=None):
    '''get velocity for specify frames with given mask

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str, default None (use all atoms), optional
        atom mask
    frame_indices : iterable that produces integer, default None, optional
        if not None, only get velocity for given frame indices

    Return
    ------
    out : 3D numpy array, shape (n_frames, n_atoms, 3)

    Examples
    --------
    >>> vels = pt.get_velocity(traj, frame_indices=[0, 3]) # doctest: +SKIP

    Notes
    -----
    Since pytraj has limited support for force and velocity info, if user wants to
    load both from disk, should iterate the TrajectoryIterator (got by pytraj.iterload method)

    .. code-block:: python

        import pytraj as pt
        forces = []
        velocities = []

        traj = pt.iterload(filename, topology_filename)

        for frame in traj:
            forces.append(frame.force)
            velocities.append(frame.velocity)

        # Note: pytraj can efficiently load arbitary frame indices to memory
        for frame in pt.iterframe(traj, frame_indices=[0, 8, 8, 100, 1000]): pass
    '''
    if mask is None:
        atm_indices = None
    else:
        if not isinstance(mask, string_types):
            # array-like
            atm_indices = mask
        else:
            atm_indices = traj.top.select(mask)

    fi = traj.iterframe(frame_indices=frame_indices)
    n_atoms = traj.n_atoms if mask is None else len(atm_indices)
    n_frames = fi.n_frames

    data = np.empty((n_frames, n_atoms, 3), dtype='f8')
    for idx, frame in enumerate(fi):
        if not frame.has_velocity():
            raise ValueError('frame does not have velocity')
        data[idx] = frame.velocity if mask is None else frame.velocity[
            atm_indices]
    return data


@super_dispatch()
def randomize_ions(traj, mask, around, by, overlap, seed=1, top=None, frame_indices=None):
    """randomize_ions for given Frame with Topology

    Parameters
    ----------
    traj : Trajectory-like or a Frame
        ``traj`` must be mutable
    mask : str
        cpptraj command
    frame_indices : {None, array-like}, optional
    top : Topology, optional (only needed if ``traj`` does not have it)

    Examples
    --------
    >>> pt.randomize_ions(traj, mask='@Na+', around=':1-16', by=5.0, overlap=3.0, seed=113698) # doctest: +SKIP
    """
    around_ = 'around ' + str(around)
    by_ = 'by ' + str(by)
    overlap_ = 'overlap ' + str(overlap)
    seed_ = 'seed ' + str(seed)
    command = ' '.join((mask, around_, by_, overlap_, seed_))
    _assert_mutable(traj)

    act = c_action.Action_RandomizeIons()
    act(command, traj, top=top)
    return traj


def clustering_dataset(array_like, command=''):
    '''cluster dataset

    Returns
    -------
    cluster index for each data point

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np
    >>> array_like = np.random.randint(0, 10, 1000)
    >>> data = pt.clustering_dataset(array_like, 'clusters 10 epsilon 3.0')
    '''
    c_dslist = CpptrajDatasetList()
    c_dslist.add('double', '__array_like')
    c_dslist[0].resize(len(array_like))
    c_dslist[0].values[:] = array_like
    act = c_analysis.Analysis_Clustering()
    command = 'data __array_like ' + command
    act(command, dslist=c_dslist)

    return np.array(c_dslist[-1])


@register_pmap
@super_dispatch()
def multidihedral(traj=None,
                  dhtypes=None,
                  resrange=None,
                  define_new_type=None,
                  range360=False,
                  dtype='dataset',
                  top=None,
                  frame_indices=None):
    """perform dihedral search

    Parameters
    ----------
    traj : Trajectory-like object
    dhtypes : dihedral type, default None
        if None, calculate all supported dihedrals
    resrange : str | array-like
        residue range for searching. If `resrange` is string, use index starting with 1
        (cpptraj convertion)
        if `resrange` is array-like, use index starting from 0 (python convention)
    define_new_type : str
        define new type for searching
    range360 : bool, default False
        if True: use 0-360
    top : Topology | str, optional
        only need to have 'top' if can not find it in `traj`


    Returns
    -------
    pytraj.DatasetList (use `values` attribute to get raw `numpy` array)

    Notes
    -----
        Dataset lables show residue number in 1-based index

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.multidihedral(traj)
    >>> data = pt.multidihedral(traj, 'phi psi')
    >>> data = pt.multidihedral(traj, resrange=range(8))
    >>> data = pt.multidihedral(traj, range360=True)
    >>> data = pt.multidihedral(traj, resrange='1,3,5')
    >>> data = pt.multidihedral(traj, dhtypes='phi psi')
    >>> data = pt.multidihedral(traj, dhtypes='phi psi', resrange='3-7')
    >>> data = pt.multidihedral(traj, dhtypes='phi psi', resrange=[3, 4, 8])
    """
    if resrange:
        if isinstance(resrange, string_types):
            _resrange = "resrange " + str(resrange)
        else:
            from pytraj.utils import convert as cv
            _resrange = cv.array_to_cpptraj_range(resrange)
            _resrange = "resrange " + str(_resrange)
    else:
        _resrange = " "

    if dhtypes:
        d_types = str(dhtypes)
    else:
        d_types = " "

    if define_new_type:
        dh_types = ' '.join(('dihtype', str(define_new_type)))
    else:
        dh_types = ''

    if range360:
        _range360 = 'range360'
    else:
        _range360 = ''

    command_ = " ".join((d_types, _resrange, dh_types, _range360))

    c_dslist = CpptrajDatasetList()
    act = c_action.Action_MultiDihedral()
    act(command_, traj, top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_multidihedral = multidihedral


@super_dispatch()
def atomicfluct(traj=None,
                mask="",
                top=None,
                dtype='ndarray',
                frame_indices=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.atomicfluct(traj, '@CA')
    >>> data[:3]
    array([[  5.        ,   0.61822273],
           [ 16.        ,   0.5627449 ],
           [ 40.        ,   0.53717119]])
    '''
    c_dslist = CpptrajDatasetList()
    act = c_action.Action_AtomicFluct()
    act(mask, traj, top=top, dslist=c_dslist)
    act.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_atomicfluct = atomicfluct


def bfactors(traj=None,
             mask="",
             byres=True,
             top=None,
             dtype='ndarray',
             frame_indices=None):
    # Not: do not use super_dispatch here since we used in calc_atomicfluct
    """calculate pseudo bfactor

    Notes
    -----
    This is **NOT** getting bfactor from xray, but computing bfactor from simulation.

    Parameters
    ----------
    traj: Trajectory-like
    mask: str, mask

    Returns
    -------
    if dtype is 'ndarray' (default), return a numpy array
    with shape=(n_atoms/n_residues, 2) ([atom_or_residue_idx, value])

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.load(fn, tn, mask='!:WAT')
    >>> traj = pt.superpose(traj)
    >>> bfactor = pt.calc_bfactors(traj, byres=True)
    """
    byres_text = "byres" if byres else ""

    # need to convert to string mask
    # do not use super_dispatch again
    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)
    command_ = " ".join((mask, byres_text, "bfactor"))
    return calc_atomicfluct(traj=traj,
                            mask=command_,
                            top=top,
                            dtype=dtype,
                            frame_indices=frame_indices)

calc_bfactors = bfactors


@register_pmap
def vector(traj=None,
           command="",
           frame_indices=None,
           dtype='ndarray',
           top=None):
    """perform vector calculation. See example below. Same as 'vector' command in cpptraj.

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces :class:`pytraj.Frame`
    command : str or a list of strings, cpptraj command
    frame_indices : array-like, optional, default None
        only perform calculation for given frame indices
    dtype : output's dtype, default 'ndarray'
    top : Topology, optional, default None

    Returns
    -------
    out : numpy ndarray, shape (n_frames, 3) if command is a string
          numpy ndarray, shape (n_vectors, n_frames, 3) if command is a list of strings

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.calc_vector(traj, "@CA @CB")
    >>> data = pt.calc_vector(traj, [("@CA @CB"),])
    >>> data = pt.calc_vector(traj, "principal z")
    >>> data = pt.calc_vector(traj, "principal x")
    >>> data = pt.calc_vector(traj, "ucellx")
    >>> data = pt.calc_vector(traj, "boxcenter")
    >>> data = pt.calc_vector(traj, "box")

    Notes
    -----
    It's faster to calculate with a list of commands.
    For example, if you need to perform 3 calculations for 'ucellx', 'boxcenter', 'box'
    like below:

    >>> data = pt.calc_vector(traj, "ucellx")
    >>> data = pt.calc_vector(traj, "boxcenter")
    >>> data = pt.calc_vector(traj, "box")

    You should use a list of commands for faster calculation.

    >>> comlist = ['ucellx', 'boxcenter', 'box']
    >>> data = pt.calc_vector(traj, comlist)
    """

    c_dslist = CpptrajDatasetList()
    top_ = get_topology(traj, top)
    list_of_commands = get_list_of_commands(command)
    fi = get_fiterator(traj, frame_indices)
    actlist = ActionList()

    for command in list_of_commands:
        act = c_action.Action_Vector()
        actlist.add(act, command, top_, dslist=c_dslist)
    actlist.compute(fi)

    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_vector = vector


@super_dispatch()
def _calc_vector_center(traj=None,
                        mask="",
                        top=None,
                        mass=False,
                        dtype='ndarray',
                        frame_indices=None):

    c_dslist = CpptrajDatasetList()
    c_dslist.set_own_memory(False)  # need this to avoid segmentation fault
    act = c_action.Action_Vector()
    command = "center " + mask

    if mass:
        command += " mass"

    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype=dtype)


@register_pmap
def center_of_mass(traj=None,
                   mask='',
                   top=None,
                   dtype='ndarray',
                   frame_indices=None):
    '''compute center of mass

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # compute center of mass residue 3 for first 2 frames.
    >>> pt.calc_center_of_mass(traj(stop=2), ':3')
    array([[-0.661702  ,  6.69124347,  3.35159413],
           [ 0.5620708 ,  7.82263042, -0.72707798]])
    '''
    # note: do not use super_dispatch for this method since
    # we already use for _calc_vector_center
    return _calc_vector_center(traj=traj,
                               mask=mask,
                               top=top,
                               mass=True,
                               dtype=dtype,
                               frame_indices=frame_indices)


calc_center_of_mass = center_of_mass


@register_pmap
@super_dispatch()
def center_of_geometry(traj=None,
                       mask="",
                       top=None,
                       dtype='ndarray',
                       frame_indices=None):

    atom_mask_obj = top(mask)
    c_dslist = CpptrajDatasetList()
    c_dslist.add("vector")

    for frame in iterframe_master(traj):
        c_dslist[0].append(frame.center_of_geometry(atom_mask_obj))
    return get_data_from_dtype(c_dslist, dtype=dtype)


calc_center_of_geometry = center_of_geometry


# do not use super_dispatch here since we did in inside this method
# to avoid complicated code checking.

@register_openmp
def pairwise_rmsd(traj=None,
                  mask="",
                  metric='rms',
                  top=None,
                  dtype='ndarray',
                  mat_type='full',
                  frame_indices=None):
    """calculate pairwise rmsd with different metrics.

    Parameters
    ----------
    traj : Trajectory-like or iterable object
    mask : mask
        if mask is "", use all atoms
    metric : {'rms', 'dme', 'srmsd', 'nofit'}
        if 'rms', perform rms fit
        if 'dme', use distance RMSD
        if 'srmsd', use symmetry-corrected RMSD
        if 'nofit', perform rmsd without fitting
    top : Topology, optional, default=None
    dtype: ndarray
        return type
    mat_type : str, {'full', 'half'}
        if 'full': return 2D array, shape=(n_frames, n_frames)
        if 'half': return 1D array, shape=(n_frames*(n_frames-1)/2, )

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> arr = pt.pairwise_rmsd(traj(0, 1000, mask='@CA'))

    >>> # calculate pairwise rmsd for all frames using CA atoms, use `dme` (distance RMSD)
    >>> # convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="dme", dtype='ndarray')

    >>> # calculate pairwise rmsd for all frames using CA atoms, nofit for RMSD
    >>> # convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="nofit", dtype='ndarray')

    >>> # calculate pairwise rmsd for all frames using CA atoms
    >>> # use symmetry-corrected RMSD, convert to numpy array
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="srmsd", dtype='ndarray')

    >>> # use different dtype
    >>> arr_np = pt.pairwise_rmsd(traj, "@CA", metric="srmsd", dtype='dataset')

    Notes
    -----
    Install ``libcpptraj`` with ``openmp`` to get benefit from parallel
    """
    # we copy Frame coordinates to DatasetCoordsCRD first

    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)

    act = c_analysis.Analysis_Rms2d()

    crdname = 'default_coords'
    c_dslist, top_, command = get_fi_with_dslist(traj, mask, frame_indices, top, crdname=crdname)

    command = ' '.join((command, metric, "crdset {} rmsout mycrazyoutput".format(crdname)))

    act(command, dslist=c_dslist)
    # remove dataset coords to free memory
    c_dslist.remove_set(c_dslist[0])

    if dtype == 'ndarray':
        return get_matrix_from_dataset(c_dslist[0], mat_type)
    else:
        return get_data_from_dtype(c_dslist, dtype)

calc_pairwise_rmsd = pairwise_rmsd


@register_pmap
def rmsd_perres(traj=None,
                mask="",
                ref=0,
                mass=False,
                resrange=None,
                perres_mask=None,
                perres_center=False,
                perres_invert=False,
                frame_indices=None,
                top=None,
                dtype='dataset', **kwd):
    """superpose ``traj`` to ``ref`` with `mask`, then calculate nofit rms for residues
    in `resrange` with given `perresmask`

    Returns
    -------
    out : pytraj.DatasetList, shape=(1+n_residues, n_frames)
        out[0]: regular rmsd
        out[1:]: perres rmsd for all given residues
        `out.values` will return corresponding numpy array
    """
    range_ = 'range %s ' % resrange
    perresmask_ = 'perresmask ' + perres_mask if perres_mask is not None else ''
    perrestcenter_ = 'perrescenter' if perres_center else ''
    perrestinvert_ = 'perresinvert' if perres_invert else ''

    cm = " ".join((mask, 'perres', range_, perresmask_, perrestcenter_,
                   perrestinvert_))
    return rmsd(traj=traj,
                mask=cm,
                ref=ref,
                nofit=False,
                mass=mass,
                frame_indices=frame_indices,
                top=top,
                dtype=dtype, **kwd)


@register_pmap
def rmsd_nofit(traj=None,
               mask="",
               ref=0,
               mass=False,
               frame_indices=None,
               top=None,
               dtype='ndarray', **kwd):
    '''compute rmsd without fitting (translating and rotating)

    Parameters
    ----------
    traj : Trajectory-like
    mask : str
    ref : Frame or int
    mass : bool, default False
        if True, use mass-weighted
    frame_indices : 1D array-like, default None
        if given, only perform calculation for those frames

    Notes
    -----

    This method is equal to pytraj.rmsd(traj, mask, ref, nofit=True, ...)

    '''
    return rmsd(traj=traj,
                mask=mask,
                ref=ref,
                mass=mass,
                nofit=True,
                frame_indices=frame_indices,
                top=top,
                dtype=dtype, **kwd)

calc_rmsd_nofit = rmsd_nofit


@register_pmap
def rmsd(traj=None,
         mask="",
         ref=0,
         ref_mask='',
         nofit=False,
         mass=False,
         update_coordinate=True,
         frame_indices=None,
         top=None,
         dtype='ndarray'):
    """compute rmsd

    Parameters
    ----------
    traj : Trajectory-like
    mask : str or 1D array-like of string or 1D or 2D array-like
        Atom mask/indices
    ref : {Frame, int}, default=0 (first frame)
        Reference frame or index.
    ref_mask: str, optional
        if given, use it instead of `mask`
    nofit : bool, default False
        if False, perform fitting (rotation and translation).
        if ``traj`` is mutable, its coordinates will be updated
        if True, not fitting.
    mass : bool, default False
        if True, include mass
    update_coordinate : bool, default True
        if True, coordinates will be updated. But this only apply to mutable Trajectory
        if False (same as `nomod` in cpptraj), no modification
    frame_indices : int 1D array-like, default None
        if not None, only compute rmsd for given frame indices
    top : {Topology, str}, default None, optional
    dtype : return data type, default='ndarray'

    Notes
    -----

    - if traj and ref has diffrent n_atoms, make sure to update ref.top
    - you can use `pytraj.rmsd` to superpose structure (use update_coordinate=True)


    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_trpcage()
    >>> # all atoms, do fitting, using ref=traj[-3]
    >>> data = pt.rmsd(traj, ref=-3)

    >>> # rmsd for 3 maskes, do fitting, using ref=traj[0] (defaul)
    >>> data = pt.rmsd(traj, mask=['@CA', '@C', ':3-18@CA'], dtype='dataset')

    >>> # rmsd to first frame, use mass ':3-13' but do not perorm fitting
    >>> data= pt.rmsd(traj, ref=traj[0], mask=':3-13', nofit=True)

    >>> # use atom indices for mask
    >>> data= pt.rmsd(traj, ref=traj[0], mask=range(40), nofit=True)

    >>> # compute rmsd (and align) with reference having different atoms
    >>> trpcage_traj = pt.datafiles.load_trpcage()[:]
    >>> tz2_traj = pt.datafiles.load_tz2()[:1]
    >>> data = pt.rmsd(trpcage_traj, mask='@1-10', ref=tz2_traj, ref_mask='@11-20')
    >>> data
    array([ 2.16203842,  2.28859396,  2.15817654, ...,  2.20767189,
            2.30087764,  1.92654945])

    Notes
    -----
    if ``traj`` is mutable and update_coordinate=True, its coordinates will be updated.

    """

    nofit_ = 'nofit' if nofit else ''
    mass_ = 'mass' if mass else ''
    nomod_ = 'nomod' if not update_coordinate else ''
    options = ' '.join((nofit_, mass_, nomod_))

    if ref_mask:
        if not mask:
            raise ValueError('mask must be provided if ref_mask is given')
        if not isinstance(ref_mask, string_types):
            ref_mask = array_to_cpptraj_atommask(ref_mask)

    if isinstance(mask, string_types):
        command = [mask, ]
    else:
        try:
            cmd = np.asarray(mask)
        except ValueError:
            raise ValueError("don't mix different types")
        dname = cmd.dtype.name
        if 'str' in dname:
            command = cmd
        elif 'int' in dname:
            if cmd.ndim == 1:
                command = [array_to_cpptraj_atommask(mask), ]
            else:
                # assume ndim==2
                command = [array_to_cpptraj_atommask(x) for x in mask]
        elif 'object' in dname:
            # different array lens or mix type
            # dangerous: assume array of two array
            command = [array_to_cpptraj_atommask(x) for x in mask]
        else:
            raise ValueError("not supported")

    top_ = get_topology(traj, top)

    ref = get_reference(traj, ref)
    fi = get_fiterator(traj, frame_indices)

    alist = ActionList()
    c_dslist = CpptrajDatasetList()

    ref_top = ref.top if ref.top else top_
    c_dslist.add('reference', name='myref')
    c_dslist[-1].top = ref_top
    c_dslist[-1].add_frame(ref)

    ref_mask = ' '.join((ref_mask, 'ref myref'))

    for cm in command:
        cm_ = ' '.join((cm, ref_mask, options))
        if 'savematrices' in cm_ and dtype not in ['dataset', 'cpptraj_dataset']:
            raise ValueError('if savematrices, dtype must be "dataset"')
        alist.add(c_action.Action_Rmsd(),
                  cm_,
                  top=top_,
                  dslist=c_dslist)

    alist.compute(fi)

    # pop Reference Dataset
    c_dslist._pop(0)

    dnew = DatasetList(c_dslist)
    return get_data_from_dtype(dnew, dtype=dtype)

# alias for `calc_rmsd`
calc_rmsd = rmsd

@super_dispatch()
def align(traj,
          mask='',
          ref=0, ref_mask='',
          mass=False,
          top=None, frame_indices=None):
    """align (superpose) trajectory to given reference

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
    ref : {int, Frame}, default 0 (first frame)
    ref_mask : str, default ''
        if not given, use traj's mask
        if given, use it
    mass : Bool, default False
        if True, mass-weighted
        if False, no mas-weighted
    frame_indices : {None, array-like}, default None
       if given, only compute RMSD for those

    Examples
    --------

    Notes
    -----
    versionadded: 1.0.6
    """
    if isinstance(traj, TrajectoryIterator):
        return traj.superpose(mask=mask, ref=ref, ref_mask=ref_mask, mass=mass)
    else:
        mask_ = mask
        refmask_ = ref_mask
        reftop = ref.top if hasattr(ref, 'top') else top
        mass_ = 'mass' if mass else ''

        refname = 'myref'
        ref_command_ = 'ref {}'.format(refname)

        command = ' '.join((ref_command_, mask_, refmask_, mass_))

        if reftop is None:
            reftop = traj.top

        c_dslist = CpptrajDatasetList()
        c_dslist.add('reference', name=refname)
        c_dslist[0].top = reftop
        c_dslist[0].add_frame(ref)

        act = c_action.Action_Align()
        act.read_input(command, top=top, dslist=c_dslist)
        act.setup(top)

        for frame in traj:
            act.compute(frame)
        act.post_process()

        # remove ref
        c_dslist._pop(0)

        return traj

@super_dispatch()
def symmrmsd(traj, mask='', ref=0, ref_mask=None,
             fit=True, remap=False,
             mass=False,
             top=None, dtype='ndarray', frame_indices=None):
    """Compute symmetry-corrected RMSD

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '' (all atoms)
    ref : {int, Frame}, default 0 (first frame)
    ref_mask : {str, None}, default None
        if None, use traj's mask
        if given, use it
    fit : Bool, default True
        if True, do fitting
        if False, nofit
    mass : Bool, default False
        if True, mass-weighted
        if False, no mas-weighted
    remap : Bool, default False
        if True, frames will be modifed for symmetry as well
    dtype : str, default 'ndarray'
        return data type
    frame_indices : {None, array-like}, default None
       if given, only compute RMSD for those

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load("TYR.nc", "TYR.parm7") # doctest: +SKIP
    >>> data = pt.symmrmsd(traj, ref=0) # doctest: +SKIP

    Notes
    -----
    versionadded: 1.0.6
    """

    mask_ = mask
    refmask_ = ref_mask if ref_mask is not None else ''
    reftop = ref.top if hasattr(ref, 'top') else top
    nofit_ = 'nofit' if not fit else ''
    mass_ = 'mass' if mass else ''
    remap_ = 'remap' if remap else ''

    refname = 'myref'
    ref_command_ = 'ref {}'.format(refname)

    command = ' '.join((mask_, refmask_, nofit_, mass_, remap_, ref_command_))

    if reftop is None:
        reftop = traj.top

    c_dslist = CpptrajDatasetList()
    c_dslist.add('reference', name=refname)
    c_dslist[0].top = reftop
    c_dslist[0].add_frame(ref)

    act = c_action.Action_SymmetricRmsd()
    act.read_input(command, top=top, dslist=c_dslist)
    act.setup(top)

    for frame in traj:
        new_frame = act.compute(frame, get_new_frame=remap)
        if remap:
            frame.xyz[:] = new_frame.xyz[:]
    act.post_process()

    # remove ref
    c_dslist._pop(0)

    return get_data_from_dtype(c_dslist, dtype=dtype)

@register_pmap
@super_dispatch()
def calc_distance_rmsd(traj=None,
                       mask='',
                       ref=0,
                       top=None,
                       dtype='ndarray',
                       frame_indices=None):
    '''compute distance rmsd between traj and reference

    Parameters
    ----------
    traj : Trajectory-like or iterator that produces Frame
    ref : {int, Frame}, default 0 (1st Frame)
    mask : str
    top : Topology or str, optional, default None
    dtype : return dtype, default 'ndarray'

    Returns
    -------
    1D ndarray if dtype is 'ndarray' (default)

    Examples
    --------
    >>> import pytraj as pt
    >>> # compute distance_rmsd to last frame
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.distance_rmsd(traj, ref=-1)

    >>> # compute distance_rmsd to first frame with mask = '@CA'
    >>> data = pt.distance_rmsd(traj, ref=0, mask='@CA')
    '''
    c_dslist = CpptrajDatasetList()
    command = mask

    act = c_action.Action_DistRmsd()
    act(command, [ref, traj], top=top, dslist=c_dslist)

    # exclude ref value
    for d in c_dslist:
        d.data = d.data[1:]
    return get_data_from_dtype(c_dslist, dtype=dtype)

# alias
distance_rmsd = calc_distance_rmsd


@super_dispatch()
def align_principal_axis(traj=None, mask="*", top=None, frame_indices=None, mass=False):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (Trajectory, Frame)
    """
    _assert_mutable(traj)

    mass_ = 'mass' if mass else ''

    act = c_action.Action_Principal()
    command = ' '.join((mask, " dorotation", mass_))
    act(command, traj, top=top)
    return traj


def principal_axes(traj=None, mask='*', dorotation=False, mass=True, top=None):
    # TODO: update doc please
    """compute principal axes

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
    mass: bool, defaul True
    if `dorotation`, the system will be aligned along principal axes (apply for mutable system)
    top : Topology, optional

    Returns
    -------
    out_0 : numpy array, shape=(n_frames, 3, 3)
    out_1: numpy array with shape=(n_frames, 3)
    """
    act = c_action.Action_Principal()
    command = mask

    _dorotation = 'dorotation' if dorotation else ''
    mass_ = 'mass' if mass else ''

    if 'name' not in command:
        command += ' name pout'

    command = ' '.join((command, _dorotation, mass_))

    top_ = get_topology(traj, top)
    c_dslist = CpptrajDatasetList()
    act(command, traj, top_, dslist=c_dslist)
    return (c_dslist[0].values, c_dslist[1].values)


def _closest_iter(act, traj):
    '''

    Parameters
    ----------
    act : Action object
    traj : Trajectory-like
    '''

    for frame in iterframe_master(traj):
        new_frame = act.compute(frame, get_new_frame=True)
        yield new_frame


@register_openmp
@super_dispatch()
def closest(traj=None,
            mask='*',
            solvent_mask=None,
            n_solvents=10,
            frame_indices=None,
            dtype='iterator',
            top=None):
    """return either a new Trajectory or a frame iterator. Keep only ``n_solvents`` closest to mask

    Parameters
    ----------
    traj : Trajectory-like | list of Trajectory-like/frames | frame iterator | chunk iterator
    mask: str, default '*' (all solute atoms)
    top : Topology-like object, default=None, optional
    dtype : {'iterator', 'trajectory'}, default 'iterator'
        if 'iterator', return a tuple of Frame iterator and new Toplogy. 'iterator' is good for streaming large trajectory.
        if 'trajectory', return a new Trajectory.

    Returns
    -------
    out : (check above)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # obtain new traj, keeping only closest 100 waters
    >>> # to residues 1 to 13 (index starts from 1) by distance to the first atom of water
    >>> t = pt.closest(traj, mask='@CA', n_solvents=10)
    """
    # check if top has solvent
    c_dslist = CpptrajDatasetList()

    command = str(n_solvents) + ' ' + mask

    act = c_action.Action_Closest()

    if solvent_mask is not None:
        top = top.copy()
        top.set_solvent(solvent_mask)

    has_solvent = False
    for mol in top.mols:
        if mol.is_solvent():
            has_solvent = True
            break
    if not has_solvent:
        raise RuntimeError("Topology does not have solvent")

    act.read_input(command, top, dslist=c_dslist)
    new_top = act.setup(top, get_new_top=True)

    fiter = _closest_iter(act, traj)

    if dtype == 'trajectory':
        return Trajectory(xyz=np.array([frame.xyz.copy() for frame in fiter]), top=new_top.copy())
    else:
        # iterator
        return (fiter, new_top.copy())


@register_pmap
@super_dispatch(refindex=3)
def native_contacts(traj=None,
                    mask="",
                    mask2="",
                    ref=0,
                    dtype='dataset',
                    distance=7.0,
                    image=True,
                    include_solvent=False,
                    byres=False,
                    frame_indices=None,
                    top=None):
    """compute native contacts

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # use 1st frame as reference, don't need specify ref Frame
    >>> data = pt.native_contacts(traj)

    >>> # explicitly specify reference, specify distance cutoff
    >>> ref = traj[3]
    >>> data = pt.native_contacts(traj, ref=ref, distance=8.0)

    >>> # use integer array for mask
    >>> data = pt.native_contacts(traj, mask=range(100), mask2=[200, 201], ref=ref, distance=8.0)
    """
    ref = get_reference(traj, ref)
    act = c_action.Action_NativeContacts()
    c_dslist = CpptrajDatasetList()

    if not isinstance(mask2, string_types):
        # [1, 3, 5] to "@1,3,5
        mask2 = array_to_cpptraj_atommask(mask2)
    command = ' '.join((mask, mask2))

    _distance = str(distance)
    _noimage = "noimage" if not image else ""
    _includesolvent = "includesolvent" if include_solvent else ""
    _byres = "byresidue" if byres else ""

    command_ = " ".join(('ref myframe', command, _distance, _noimage,
                         _includesolvent, _byres))
    c_dslist.add('ref_frame', 'myframe')
    c_dslist[0].top = top
    c_dslist[0].add_frame(ref)
    act(command_, traj, top=top, dslist=c_dslist)
    c_dslist._pop(0)

    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def grid(traj=None, command="", top=None, dtype='dataset'):
    """
    """
    # TODO: doc, rename method, move to seperate module?
    act = c_action.Action_Grid()
    c_dslist = CpptrajDatasetList()

    # cpptraj require output
    command = "tmp_pytraj_grid_output.txt " + command
    with tempfolder():
        act(command, traj, dslist=c_dslist, top=top)
    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_grid = grid


@super_dispatch()
def check_structure(traj, mask='', options='', frame_indices=None, top=None, dtype='ndarray'):
    """check if the structure is ok or not

    Parameters
    ----------
    traj : Trajectory-like
    mask: str, default all atoms
    options : str, default ''
        extra cpptraj options
    dtype : str, default 'ndarray'

    Returns
    -------
    out : 1D-array
        number of problems for each frame

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()
    >>> failures = pt.check_structure(traj[0], top=traj.top)
    """
    act = c_action.Action_CheckStructure()
    command = ' '.join((mask, options))
    c_dslist = CpptrajDatasetList()

    act(command, traj, top=top, dslist=c_dslist)
    return get_data_from_dtype(c_dslist, dtype=dtype)


def timecorr(vec0,
             vec1,
             order=2,
             tstep=1.,
             tcorr=10000.,
             norm=False,
             dtype='ndarray'):
    """compute time correlation.

    Parameters
    ----------
    vec0 : 2D array-like, shape=(n_frames, 3)
    vec1 : 2D array-like, shape=(n_frames, 3)
    order : int, default 2
    tstep : float, default 1.
    tcorr : float, default 10000.
    norm : bool, default False
    """
    # TODO: doc. not yet assert to cpptraj's output
    act = c_analysis.Analysis_Timecorr()

    c_dslist = CpptrajDatasetList()

    c_dslist.add("vector", "_vec0")
    c_dslist.add("vector", "_vec1")
    c_dslist[0].data = np.asarray(vec0).astype('f8')
    c_dslist[1].data = np.asarray(vec1).astype('f8')

    order_ = "order " + str(order)
    tstep_ = "tstep " + str(tstep)
    tcorr_ = "tcorr " + str(tcorr)
    norm_ = "norm" if norm else ""
    command = " ".join(('vec1 _vec0 vec2 _vec1', order_, tstep_, tcorr_, norm_))
    act(command, dslist=c_dslist)
    return get_data_from_dtype(c_dslist[2:], dtype=dtype)

@super_dispatch()
def velocityautocorr(traj, mask='', maxlag=-1, tstep=1.0, direct=True, norm=False,
        usevelocity=False,
        dtype='ndarray',
        top=None,
        velocity_arr=None,
        ):
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
    usevelocity : bool, default False
        if True, use velocity info in Frame
    dtype : str, default 'ndarray'
        return data type
    top : None or Topology, default None, optional
    velocity_arr : None or 3D like array, default None
        only use `velocity_arr` if usevelocity is True

    Notes
    -----
    If you create Trajectory by `pytraj.load` method, there is no velocity information.
    So if you want to use `usevelocity=True`, you need to provide 3D-array velocity_arr
    """
    from pytraj import Frame

    if velocity_arr is not None:
        velocity_arr = np.asarray(velocity_arr)

        if len(velocity_arr.shape) != 3:
            raise ValueError('provided velocity_arr must be 3D array-like, shape=(n_frames, n_atoms, 3)')

    act = c_action.Action_VelocityAutoCorr()
    c_dslist = CpptrajDatasetList()

    maxlag_ = ' '.join(('maxlag', str(maxlag)))
    tstep_ = ' '.join(('tstep', str(tstep)))
    direct_ = 'direct' if direct else ''
    norm_ = 'norm' if norm else ''
    usevelocity_ = 'usevelocity' if usevelocity else ''

    command = ' '.join((maxlag_, tstep_, direct_, norm_, usevelocity_))
    crdinfo = dict(has_velocity=True) if usevelocity else dict()

    act.read_input(command, top, dslist=c_dslist)
    act.setup(top, crdinfo=crdinfo)

    frame_template = Frame()

    if usevelocity and velocity_arr is not None:
        frame_template._allocate_force_and_velocity(top, crdinfo=dict(has_velocity=True))
        use_template = True
    else:
        use_template = False

    for idx, frame in enumerate(traj):
        if not use_template:
            if usevelocity and not frame.has_velocity():
                raise ValueError("Frame must have velocity if specify 'usevelocity'")
            act.compute(frame)
        else:
            vel = velocity_arr[idx]
            frame_template.xyz[:] = frame.xyz[:]
            frame_template.velocity[:] = vel
            act.compute(frame_template)

    act.post_process()

    return get_data_from_dtype(c_dslist, dtype=dtype)


def crank(data0, data1, mode='distance', dtype='ndarray'):
    """

    Parameters
    ----------
    data0 : array-like
    data1 : array-like
    mode : str, {'distance', 'angle'}
    dtype : str

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> distances = pt.distance(traj, [':3 :7', ':8 :12'])
    >>> out = pt.crank(distances[0], distances[1])

    Notes
    -----
    Same as `crank` in cpptraj
    """
    c_dslist = CpptrajDatasetList()
    c_dslist.add("double", "d0")
    c_dslist.add("double", "d1")

    c_dslist[0].data = np.asarray(data0)
    c_dslist[1].data = np.asarray(data1)

    act = c_analysis.Analysis_CrankShaft()
    command = ' '.join((mode, 'd0', 'd1'))
    act(command, dslist=c_dslist)
    return get_data_from_dtype(c_dslist[2:], dtype=dtype)


@super_dispatch()
def search_neighbors(traj=None,
                     mask='',
                     frame_indices=None,
                     dtype='dataset',
                     top=None):
    """search neighbors

    Returns
    -------
    :ref:`pytraj.DatasetList`, is a list of atom index arrays for each frame.
    Those arrays might not have the same lenghth

    Note
    ----
    not validate yet

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> indices = pt.search_neighbors(traj, ':5<@5.0') # around residue 5 with 5.0 cutoff
    """
    c_dslist = DatasetList()

    for idx, frame in enumerate(iterframe_master(traj)):
        top.set_distance_mask_reference(frame)
        c_dslist.append({str(idx): np.asarray(top.select(mask))})
    return get_data_from_dtype(c_dslist, dtype)


@register_pmap
def pucker(traj=None,
           pucker_mask=("C1'", "C2'", "C3'", "C4'", "O4'"),
           resrange=None,
           top=None,
           dtype='dataset',
           range360=False,
           method='altona',
           use_com=True,
           amplitude=False,
           offset=None):
    """compute pucker

    Parameters
    ----------
    traj : Trajectory-like
    pucker_mask : str
    resrange : None or array of int
    top : Topology, optional
    dtype : str, return type
    range360: bool, use 360 or 180 scale
    method : {'altona', 'cremer'}, default 'altona'
    use_com : bool
    amplitude : bool, default False
    offset : None or float

    Returns
    -------
    Dataset
    """
    top_ = get_topology(traj, top)
    if resrange is None:
        resrange = range(top_.n_residues)

    _range360 = "range360" if range360 else ""
    geom = "geom" if not use_com else ""
    amp = "amplitude" if amplitude else ""
    offset_ = "offset " + str(offset) if offset else ""

    c_dslist = CpptrajDatasetList()

    for res in resrange:
        command = " ".join((":" + str(res + 1) + '@' + x for x in pucker_mask))
        name = "pucker_res" + str(res + 1)
        command = " ".join((name, command, _range360, method, geom, amp,
                            offset_))

        act = c_action.Action_Pucker()
        act(command, traj, top=top_, dslist=c_dslist)

    return get_data_from_dtype(c_dslist, dtype)

calc_pucker = pucker


@super_dispatch()
def center(traj=None,
           mask="",
           center='box',
           mass=False,
           top=None,
           frame_indices=None):
    """center

    Parameters
    ----------
    traj : Trajectory-like or Frame iterator
    mask : str, mask
    center : str, {'box', 'origin'}
    mass : bool, default: False
        if True, use mass weighted
    top : Topology, optional, default: None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # load all frames to memory so we can 'mutate' them
    >>> traj = traj[:]
    >>> # all atoms, center to box center (x/2, y/2, z/2)
    >>> traj = pt.center(traj)

    >>> # center at origin, use @CA
    >>> traj = pt.center(traj, '@CA', center='origin')

    >>> # center to box center, use mass weighted
    >>> traj = pt.center(traj, mass=True)
    >>> traj = pt.center(traj, ':1', mass=True)

    Returns
    -------
    updated traj

    See also
    --------
    pytraj.translate
    """
    if center.lower() not in ['box', 'origin']:
        raise ValueError('center must be box or origin')
    center_ = '' if center == 'box' else center
    mass_ = 'mass' if mass else ''
    command = ' '.join((mask, center_, mass_))
    _assert_mutable(traj)

    act = c_action.Action_Center()
    act(command, traj, top=top)
    return traj


def rotate_dihedral(traj=None, mask="", top=None):
    # change to pt.rotate_dihedral(traj, res=0,
    #              mask=("O4'", "C1'", "N9", "C4"), deg=120)?
    """
    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()[:]
    >>> traj = pt.rotate_dihedral(traj, "3:chin:120") # rotate chin of res 3 to 120 deg
    >>> traj = pt.rotate_dihedral(traj, "1:O4':C1':N9:C4:120") # rotate dihedral with given mask

    Returns
    -------
    updated traj

    Notes
    -----
    Syntax and method's name might be changed
    """
    _assert_mutable(traj)
    top_ = get_topology(traj, top)

    if "custom:" in mask:
        command = mask
    else:
        command = "custom:" + mask

    act = c_action.Action_MakeStructure()

    act(command, traj, top=top_)
    return traj


@register_openmp
def replicate_cell(traj=None, mask="", direction='all', top=None):
    '''create a trajectory where the unit cell is replicated in 1 or more direction (up to 27)

    Parameters
    ----------
    traj : Trajectory-like or Frame iterator
    mask : str, default: ""
        if default, using all atoms
        else: given mask
    direction: {'all', 'dir'} or list/tuple of <XYZ> (below)
        if 'all', replicate cell once in all possible directions
        if 'dir', need to specify the direction with format 'dir <XYZ>', where each X (Y, Z)
        is either 0, 1 or -1 (see example below)
    top : Topology, optional, default: None

    Returns
    -------
    traj : pytraj.Trajectory

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> new_traj = pt.replicate_cell(traj, direction='all')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 111')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 1-10')
    >>> new_traj = pt.replicate_cell(traj, direction='dir 001 dir 1-10')

    >>> # similiar usage
    >>> new_traj = pt.replicate_cell(traj, direction=('001', '0-10'))
    >>> new_traj = pt.replicate_cell(traj, direction=['001', '0-10'])
    '''
    top_ = get_topology(traj, top)
    if isinstance(direction, string_types):
        _direction = direction
    elif isinstance(direction, (list, tuple)):
        # example: direction = ('001, '0-10')
        _direction = 'dir ' + ' dir '.join(direction)
    else:
        raise ValueError(
            'only support ``direction`` as a string or list/tuple of strings')
    command = ' '.join(('name tmp_cell', _direction, mask))

    act = c_action.Action_ReplicateCell()
    c_dslist = CpptrajDatasetList()
    act(command, traj, top=top_, dslist=c_dslist)
    traj = Trajectory(xyz=c_dslist[0].xyz, top=c_dslist[0].top)

    return traj


def set_dihedral(traj, resid='1', dihedral_type=None, deg=0, top=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # make mutable traj by loading all frames to disk
    >>> traj = traj[:]
    >>> traj = pt.set_dihedral(traj, resid='3', dihedral_type='phi', deg=60)
    >>> traj = pt.set_dihedral(traj, resid=2, dihedral_type='phi', deg=60)

    Returns
    -------
    updated traj
    '''
    top_ = get_topology(traj, top)

    if not isinstance(resid, string_types):
        resid = str(resid + 1)
    deg = str(deg)

    command = ':'.join((dihedral_type, resid, dihedral_type, deg))
    make_structure(traj, command, top=top_)
    return traj


def make_structure(traj=None, mask="", top=None):
    _assert_mutable(traj)
    top_ = get_topology(traj, top)

    command = mask
    act = c_action.Action_MakeStructure()
    act(command, traj, top=top_)


@super_dispatch()
def projection(traj,
               mask='',
               eigenvalues=None,
               eigenvectors=None,
               scalar_type='covar',
               average_coords=None,
               frame_indices=None,
               dtype='ndarray',
               top=None):
    '''compute projection along given eigenvectors

    Parameters
    ----------
    traj : Trajectory-like
    mask : atom mask, either string or array-like
    eigenvalues : 1D array-like
    eigenvectors : 2D array-like
    scalar_type : str, {'covar', 'mwcovar', }, default 'covar'
        make sure to provide correct scalar_type.
        Note: not yet support 'dihcovar' and 'idea'
    average_coords : 3D array-like, optional, default None
        average coordinates. If 'None', pytraj will compute mean_structure with given mask
    frame_indices : array-like
        If not None, compute projection for given frames.
    dtype : str, return data type, default 'ndarray'
    top : Topology, optional, default None

    Returns
    -------
    projection_data : ndarray, shape=(n_vecs, n_frames)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> mat = pt.matrix.covar(traj, '@CA')
    >>> eigenvalues, eigenvectors = pt.matrix.diagonalize(mat, 2)

    >>> # since we compute covariance matrix, we need to specify
    >>> # scalar_type = 'covar'
    >>> scalar_type = 'covar'
    >>> data = pt.projection(traj, '@CA', eigenvalues=eigenvalues, eigenvectors=eigenvectors, scalar_type=scalar_type)
    >>> data.shape
    (2, 101)
    '''

    act = c_action.Action_Projection()
    c_dslist = CpptrajDatasetList()

    mode_name = 'my_modes'
    c_dslist.add('modes', mode_name)

    is_reduced = False
    dataset_mode = c_dslist[-1]
    n_vectors = len(eigenvalues)
    dataset_mode._set_modes(is_reduced, n_vectors, eigenvectors.shape[1],
                            eigenvalues, eigenvectors.flatten())
    dataset_mode.scalar_type = scalar_type

    if average_coords is not None:
        dataset_mode._allocate_avgcoords(3 * average_coords.shape[0])
        dataset_mode._set_avg_frame(average_coords.flatten())
    else:
        frame = mean_structure(traj, mask)
        average_coords = frame.xyz
        dataset_mode._allocate_avgcoords(3 * average_coords.shape[0])
        dataset_mode._set_avg_frame(average_coords.flatten())

    mask_ = mask
    evecs_ = 'evecs {}'.format(mode_name)
    beg_end_ = 'beg 1 end {}'.format(n_vectors)

    command = ' '.join((evecs_, mask_, beg_end_))
    act(command, traj, top=top, dslist=c_dslist)

    c_dslist._pop(0)

    return get_data_from_dtype(c_dslist, dtype=dtype)


def pca(traj,
        mask,
        n_vecs=2,
        fit=True,
        ref=None,
        ref_mask=None,
        dtype='ndarray',
        top=None):
    '''perform PCA analysis by following below steps:

    - (optional) perform rmsfit to reference if needed
    - compute covariance matrix
    - diagonalize the matrix to get eigenvectors and eigenvalues
    - perform projection of each frame with mask to each eigenvector

    Parameters
    ----------
    traj : Trajectory
        traj must be ``pytraj.Trajectory``, which can be created by ``pytraj.load`` method.
    mask : str
        atom mask for covariance matrix and projection
    n_vecs : int, default 2
        number of eigenvectors. If user want to compute projection for all eigenvectors,
        should specify n_vecs=-1 (or a negative number)
    fit : bool, default True
        if True, perform fitting before compute covariance matrix
        if False, no fitting (keep rotation and translation). In this case, `pytraj` will ignore `ref` argument.
    ref : {None, Frame, int}, default None
        if None, trajectory will be superposed to average structure
        if is Frame or integer value, trajectory will be superposed to given reference
    ref_mask : {None, str}, default None (use `mask`)
        if None, use `mask` for fitting
        if str, use this given mask for fitting
    dtype : return datatype
    top : Topology, optional

    Returns
    -------
    out1: projection_data, ndarray with shape=(n_vecs, n_frames)
    out2: tuple of (eigenvalues, eigenvectors)

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()[:]

    >>> # compute pca for first and second modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=2)
    >>> # get projection values
    >>> pca_data[0] # doctest: +SKIP
    array([[  4.93425131,  13.80002308,  20.61605835, ..., -57.92280579,
            -61.25728607, -52.85142136],
           [  4.03333616,  -6.9132452 , -14.53991318, ...,  -6.757936  ,
              2.1086719 ,  -3.60922861]], dtype=float32)
    >>> # get eigenvalues for first 2 modes
    >>> pca_data[1][0] # doctest: +SKIP
    array([ 1399.36472919,   240.42342439])

    >>> # compute pca for all modes
    >>> pca_data = pt.pca(traj, '!@H=', n_vecs=-1)

    >>> # does not perform fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=False)

    >>> # provide different mask for fitting
    >>> data = pt.pca(traj, mask='!@H=', fit=True, ref=0, ref_mask='@CA')
    '''
    # TODO: move to another file
    # NOTE: do not need to use super_dispatch here since we already use in projection
    from pytraj import matrix

    ref_mask_ = ref_mask if ref_mask is not None else mask

    if not isinstance(traj, (Trajectory, TrajectoryIterator)):
        raise ValueError('must be Trajectory-like')

    if fit:
        if ref is None:
            traj.superpose(ref=0, mask=ref_mask_)
            avg = mean_structure(traj)
            traj.superpose(ref=avg, mask=ref_mask_)
            n_refs = 2
        else:
            ref = get_reference(traj, ref)
            traj.superpose(ref=ref, mask=ref_mask_)
            n_refs = 1

    avg2 = mean_structure(traj, mask=mask)

    mat = matrix.covar(traj, mask)
    if n_vecs < 0:
        n_vecs = mat.shape[0]
    else:
        n_vecs = n_vecs

    eigenvalues, eigenvectors = matrix.diagonalize(mat, n_vecs=n_vecs, dtype='tuple')
    projection_data = projection(traj, mask=mask, average_coords=avg2.xyz,
                                 eigenvalues=eigenvalues,
                                 eigenvectors=eigenvectors,
                                 scalar_type='covar', dtype=dtype)

    # release added transformed commands for TrajectoryIterator
    if fit and hasattr(traj, '_transform_commands'):
        for _ in range(n_refs):
            traj._transform_commands.pop()
        if traj._transform_commands:
            traj._reset_transformation()
        else:
            traj._remove_transformations()
    return projection_data, (eigenvalues, eigenvectors)

calc_pca = pca


@register_openmp
@super_dispatch()
def atomiccorr(traj,
               mask='',
               cut=None,
               min_spacing=None,
               byres=True,
               frame_indices=None,
               dtype='ndarray',
               top=None):
    '''compute average correlations between the motion of atoms in mask.

    Parameters
    ----------
    traj : Trajectory-like
    mask : atom mask
    cut : {None, float}, default None
        if not None, only print correlations with absolute value greater than cut
    min_spacing : {None, float}, default None
        if not None, only calculate correlations for motion vectors spaced min_spacing apart
    byres : bool, default False
        if False, compute atomic motion vetor
        if True, Calculate motion vectors for entire residues (selected atoms in residues only).
    '''
    _mask = 'out tmp.dat ' + mask
    _cut = 'cut ' + str(cut) if cut is not None else ''
    _min_spacing = 'min ' + str(min_spacing) if min_spacing is not None else ''
    _byres = 'byres' if byres else 'byatom'
    command = ' '.join((_mask, _cut, _min_spacing, _byres))

    act = c_action.Action_AtomicCorr()
    c_dslist = CpptrajDatasetList()

    with tempfolder():
        act(command, traj, top=top, dslist=c_dslist)
        # need to post_process for this Action
        act.post_process()

    return get_data_from_dtype(c_dslist, dtype=dtype)

calc_atomiccorr = atomiccorr


def _grid(traj,
          mask,
          grid_spacing,
          offset=1,
          frame_indices=None,
          dtype='ndarray',
          top=None):
    # TODO: what's about calc_grid?
    '''make grid for atom in mask

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, atom mask
    grid_spacing : array-like, shape=(3,)
        grid spacing in X/Y/Z directions
    offset : int, optional
        bin offset, number of bins to add to each direction to grid
    dtype : str, default 'ndarray'
        output data type
    '''
    top_ = get_topology(traj, top)
    fi = get_fiterator(traj, frame_indices)
    act = c_action.Action_Bounds()
    c_dslist = CpptrajDatasetList()
    dx, dy, dz = grid_spacing
    dx_ = 'dx ' + str(dx) if dx > 0. else ''
    dy_ = 'dy ' + str(dy) if dy > 0. else ''
    dz_ = 'dz ' + str(dz) if dz > 0. else ''
    offset_ = 'offset ' + str(offset)
    command = ' '.join((mask, 'out tmp_bounds.dat', dx_, dy_, dz_,
                        'name grid_', offset_))

    with tempfolder():
        act(command, fi, top=top_, dslist=c_dslist)
    act.post_process()

    return get_data_from_dtype(c_dslist, dtype=dtype)


def transform(traj, by, frame_indices=None):
    '''transform pytraj.Trajectory by a series of cpptraj's commands

    Parameters
    ----------
    traj : Mutable Trajectory
    by : list of cpptraj commands
    frame_indices : {None, array-like}, default None
        if not None, perform tranformation for specific frames.

    Returns
    -------
    transformed Trajectory. Trajectory's coordinates will be inplace-updated

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> # perform 'autoimage', then 'rms', then 'center'
    >>> traj = pt.transform(traj[:], by=['autoimage', 'rms', 'center :1-5'])
    '''
    return traj.transform(by, frame_indices=frame_indices)


def lowestcurve(data, points=10, step=0.2):
    '''compute lowest curve for data

    Parameters
    ----------
    data : 2D array-like
    points : number of lowest points in each bin, default 10
    step : step size, default 0.2

    Return
    ------
    2d array
    '''
    points_ = 'points ' + str(points)
    step_ = 'step ' + str(step)
    label = 'mydata'
    command = ' '.join((label, points_, step_))

    data = np.asarray(data)

    c_dslist = CpptrajDatasetList()
    c_dslist.add('xymesh', label)
    c_dslist[0]._append_from_array(data.T)

    act = c_analysis.Analysis_LowestCurve()

    act(command, dslist=c_dslist)
    return np.array([c_dslist[-1]._xcrd(), np.array(c_dslist[-1].values)])


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
    c_dslist = CpptrajDatasetList()
    c_dslist.add("double", "d0")

    c_dslist[0].data = np.asarray(data)

    act = c_analysis.Analysis_AutoCorr()
    command = "d0 out _tmp.out"
    act(command, dslist=c_dslist)
    return get_data_from_dtype(c_dslist[1:], dtype=dtype)

auto_correlation_function = acorr


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

    c_dslist = CpptrajDatasetList()
    c_dslist.add("double", "d0")
    c_dslist.add("double", "d1")

    c_dslist[0].data = np.asarray(data0)
    c_dslist[1].data = np.asarray(data1)

    act = c_analysis.Analysis_Corr()
    act("d0 d1 out _tmp.out", dslist=c_dslist)
    return get_data_from_dtype(c_dslist[2:], dtype=dtype)

cross_correlation_function = xcorr


def superpose(traj, *args, **kwd):
    return align(traj, *args, **kwd)

def strip(obj, mask):
    '''return a new Trajectory or FrameIterator or Topology with given mask.

    Notes
    -----
    This method is trying to be smart. If you give it an in-memory Trajectory, it will
    return a corresponding in-memory Trajectory. If you give it an out-of-memory TrajectoryIterator,
    it will give you a corresponding FrameIterator object (out-of-memory). If you give it a Topology, it will
    return a new stripped Topology.

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> traj.n_atoms
    5293
    >>> pt.strip(traj, '!@CA').n_atoms
    12
    >>> pt.strip(traj.top, '!@CA').n_atoms
    12
    >>> pt.strip('!@CA', traj.top).n_atoms
    12
    >>> t0 = traj[:3]
    >>> pt.strip(t0, '!@CA').n_atoms
    12
    >>> fi = pt.iterframe(traj, stop=3)
    >>> fi = pt.strip(fi, '!@CA')
    >>> for frame in fi:
    ...    print(frame)
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

    >>> # raise ValueError
    >>> pt.strip(0, '@CA')
    Traceback (most recent call last):
        ...
    ValueError: object must be either Trajectory or Topology
    '''

    if isinstance(obj, string_types) and not isinstance(
            mask, string_types):
        obj, mask = mask, obj

    kept_mask = '!(' + mask + ')'

    if isinstance(obj, (Topology, Trajectory)):
        # return new Topology or new Trajectory
        return obj[kept_mask]
    elif isinstance(obj, TrajectoryIterator):
        # return a FrameIterator
        return obj(mask=kept_mask)
    elif hasattr(obj, 'mask'):
        obj.mask = kept_mask
        return obj
    else:
        raise ValueError('object must be either Trajectory or Topology')


def _rotdif(matrices, command):
    """

    Parameters
    ----------
    matrices : 3D array, shape=(n_frames, 3, 3)
        rotation matrices
    command : str
        full cpptraj's command

    Notes
    -----
    This method interface will be changed.
    """
    matrices = np.asarray(matrices)

    c_dslist = CpptrajDatasetList()
    c_dslist.add(dtype='matrix3x3', name='myR0')
    c_dslist[-1].aspect = "RM"
    c_dslist[-1]._append_from_array(matrices)

    command = 'rmatrix myR0[RM] ' + command
    act = c_analysis.Analysis_Rotdif()
    act(command, dslist=c_dslist)
    return get_data_from_dtype(c_dslist[1:])

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

    c_dslist = CpptrajDatasetList()
    crdname = '_DEFAULTCRD_'
    c_dslist.add('coords', name=crdname)
    c_dslist[0].top = traj.top

    for frame in traj:
        c_dslist[0].append(frame)

    act = c_analysis.Analysis_Wavelet()
    act(command, dslist=c_dslist)
    c_dslist.remove_set(c_dslist[crdname])
    return get_data_from_dtype(c_dslist, dtype='dict')
