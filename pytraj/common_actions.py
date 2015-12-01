from __future__ import absolute_import
import numpy as np

from pytraj.trajectory import Trajectory
from pytraj.trajectory_iterator import TrajectoryIterator
from ._get_common_objects import _get_topology, _get_data_from_dtype, _get_list_of_commands
from ._get_common_objects import _get_matrix_from_dataset
from ._get_common_objects import _get_reference_from_traj, _get_fiterator
from ._get_common_objects import _super_dispatch, _get_fi_with_dslist
from .utils import ensure_not_none_or_string
from .utils import is_int
from .utils.context import goto_temp_folder
from .utils.convert import array_to_cpptraj_atommask
from .externals.six import string_types
from .datasets.DatasetList import DatasetList as CpptrajDatasetList
from .datasetlist import DatasetList
from ._shared_methods import iterframe_master
from .decorators import _register_pmap, _register_openmp
from .actions import CpptrajActions
from .analyses import CpptrajAnalyses
from .core.action_list import ActionList
from .utils.convert import array2d_to_cpptraj_maskgroup

list_of_cal = ['calc_distance',
               'calc_dihedral',
               'calc_radgyr',
               'calc_angle',
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
               'calc_COM',
               'calc_center_of_mass',
               'calc_center_of_geometry',
               'calc_pairwise_rmsd',
               'calc_grid',
               'calc_temperatures',
               'calc_atomiccorr',
               'calc_bfactors',
               'calc_diffusion',
               'calc_distance_rmsd',
               'calc_mindist',
               'calc_pairwise_distance',
               'calc_rmsd_nofit',
               'calc_rotation_matrix', 
               'calc_pca',]

list_of_do = ['do_translation', 'do_rotation', 'do_autoimage', 'do_scaling']

list_of_get = ['get_average_frame', 'get_velocity']

list_of_the_rest = ['rmsd', 'align_principal_axis', 'principal_axes', 'closest',
                    'transform', 'native_contacts', 'set_dihedral',
                    'auto_correlation_function', 'cross_correlation_function',
                    'check_structure', 'mean_structure', 'lifetime', 'lowestcurve',
                    'make_structure', 'replicate_cell', 'pucker', 'rmsd_perres',
                    'timecorr', 'search_neighbors', ]

__all__ = list_of_do + list_of_cal + list_of_get + list_of_the_rest


def _2darray_to_atommask_groups(seq):
    '''

    >>> list(_2darray_to_atommask_groups([[0, 3], [4, 7]]))
    ['@1 @4', '@5 @8']
    '''
    for arr in seq:
        # example: arr = [0, 3]; turns ot '@1 @4'
        yield '@' + str(arr[0] + 1) + ' @' + str(arr[1] + 1)


def _assert_mutable(trajiter):
    if isinstance(trajiter, TrajectoryIterator):
        raise ValueError(
            "This analysis does not support immutable object. Use `pytraj.Trajectory`")


@_register_pmap
def calc_distance(traj=None,
                  mask="",
                  frame_indices=None,
                  dtype='ndarray',
                  top=None,
                  image=True,
                  n_frames=None):
    # TODO: add image, noe, ...
    """calculate distance between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or a list of string or a 2D array-like of integers
    frame_indices : array-like, optional, default None
    dtype : return type, default 'ndarray'
    top : Topology, optional
    image : bool, default True
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

    traj = _get_fiterator(traj, frame_indices)
    _top = _get_topology(traj, top)
    _noimage = 'noimage' if not image else ''

    cm_arr = np.asarray(command)

    if 'int' in cm_arr.dtype.name:

        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(cm_arr)

        if int_2darr.shape[1] != 2:
            raise ValueError("require int-array with shape=(n_atoms, 2)")

        if n_frames is None:
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")

        arr = np.empty([n_frames, len(int_2darr)])

        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._distance(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = DatasetList({'distance': arr})
            return _get_data_from_dtype(py_dslist, dtype)

    elif isinstance(command, (list, tuple, string_types, np.ndarray)):
        # create a list
        if not isinstance(command, np.ndarray):
            list_of_commands = _get_list_of_commands(command)
        else:
            list_of_commands = command

        dslist = CpptrajDatasetList()
        actlist = ActionList()

        for cm in list_of_commands:
            if not image:
                cm = ' '.join((cm, _noimage))
            actlist.add_action(CpptrajActions.Action_Distance(),
                               cm,
                               _top,
                               dslist=dslist)

        actlist.do_actions(traj)
        return _get_data_from_dtype(dslist, dtype)

    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


def calc_pairwise_distance(traj=None,
                           mask_1='',
                           mask_2='',
                           top=None,
                           dtype='ndarray',
                           frame_indices=None):
    '''calculate pairwise distance between atoms in mask_1 and atoms in mask_2

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

    _top = _get_topology(traj, top)
    indices_1 = _top.select(mask_1) if isinstance(mask_1,
                                                  string_types) else mask_1
    indices_2 = _top.select(mask_2) if isinstance(mask_2,
                                                  string_types) else mask_2
    arr = np.array(list(product(indices_1, indices_2)))
    mat = calc_distance(traj,
                        mask=arr,
                        dtype=dtype,
                        top=_top,
                        frame_indices=frame_indices)
    mat = mat.T
    return (mat.reshape(mat.shape[0], len(indices_1), len(indices_2)),
            arr.reshape(
                len(indices_1), len(indices_2), 2))


@_register_pmap
def calc_angle(traj=None,
               mask="",
               top=None,
               dtype='ndarray',
               frame_indices=None,
               *args,
               **kwd):
    """calculate angle between two maskes

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
    dslist = CpptrajDatasetList()
    act = CpptrajActions.Action_Angle()

    command = mask

    ensure_not_none_or_string(traj)

    traj = _get_fiterator(traj, frame_indices)
    _top = _get_topology(traj, top)
    cm_arr = np.asarray(command)

    if 'int' not in cm_arr.dtype.name:
        if isinstance(command, string_types):
            # need to remove 'n_frames' keyword since Action._master does not use
            # it
            try:
                del kwd['n_frames']
            except KeyError:
                pass
            # cpptraj mask for action
            act(command, traj, top=_top, dslist=dslist, *args, **kwd)
            return _get_data_from_dtype(dslist, dtype)
        elif isinstance(command, (list, tuple, np.ndarray)):
            list_of_commands = command
            dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add_action(CpptrajActions.Action_Angle(),
                                   cm,
                                   _top,
                                   dslist=dslist,
                                   *args,
                                   **kwd)
            actlist.do_actions(traj)
            return _get_data_from_dtype(dslist, dtype)

    else:
        # ndarray of integer
        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(int_2darr)

        if int_2darr.shape[1] != 3:
            raise ValueError("require int-array with shape=(n_atoms, 3)")

        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._angle(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = DatasetList({'angle': arr})
            return _get_data_from_dtype(py_dslist, dtype)


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    '''calculate dihedral within a single residue. For internal use only.

    Parameters
    ----------
    traj : Trajectory-like
    mask : tuple of strings
    resid : str, resid
    dtype

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.common_actions import _dihedral_res
    >>> traj = pt.datafiles.load_tz2()
    >>> data = _dihedral_res(traj, mask=('N', 'CA', 'C', 'O'), resid=0)
    '''

    if is_int(resid):
        resid = str(resid + 1)
    else:
        resid = resid
    m = ' :%s@' % resid
    command = m + m.join(mask)
    return calc_dihedral(traj=traj, mask=command, top=top, dtype=dtype)


@_register_pmap
def calc_dihedral(traj=None,
                  mask="",
                  top=None,
                  dtype='ndarray',
                  frame_indices=None,
                  *args,
                  **kwd):
    """calculate dihedral angle between two maskes

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
    act = CpptrajActions.Action_Dihedral()
    dslist = CpptrajDatasetList()

    ensure_not_none_or_string(traj)
    command = mask

    traj = _get_fiterator(traj, frame_indices)
    _top = _get_topology(traj, top)
    cm_arr = np.asarray(command)

    if 'int' not in cm_arr.dtype.name:
        if isinstance(command, string_types):
            # need to remove 'n_frames' keyword since Action._master does not use
            # it
            try:
                del kwd['n_frames']
            except:
                pass
            # cpptraj mask for action
            act(command, traj, top=_top, dslist=dslist)
            return _get_data_from_dtype(dslist, dtype)

        elif isinstance(command, (list, tuple, np.ndarray)):
            list_of_commands = command
            from pytraj.core.action_list import ActionList
            from pytraj.actions.CpptrajActions import Action_Dihedral
            dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add_action(Action_Dihedral(),
                                   cm,
                                   _top,
                                   dslist=dslist,
                                   *args,
                                   **kwd)
            actlist.do_actions(traj)
            return _get_data_from_dtype(dslist, dtype)
    else:
        # ndarray
        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(int_2darr)

        if int_2darr.shape[1] != 4:
            raise ValueError("require int-array with shape=(n_atoms, 4)")

        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']
        arr = np.empty([n_frames, len(int_2darr)])
        for idx, frame in enumerate(iterframe_master(traj)):
            arr[idx] = frame._dihedral(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = DatasetList({'dihedral': arr})
            return _get_data_from_dtype(py_dslist, dtype)


@_register_pmap
def calc_mindist(traj=None,
                 command="",
                 top=None,
                 dtype='ndarray',
                 frame_indices=None):
    '''calculate mindist

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.mindist(traj, '@CA @H')
    '''

    traj = _get_fiterator(traj, frame_indices)
    act = CpptrajActions.Action_NativeContacts()
    dslist = CpptrajDatasetList()

    if not isinstance(command, string_types):
        command = array2d_to_cpptraj_maskgroup(command)
    _command = "mindist " + command
    traj = _get_fiterator(traj, frame_indices)
    _top = _get_topology(traj, top)
    act(_command, traj, top=_top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)[-1]


@_super_dispatch()
def calc_diffusion(traj,
                   mask="",
                   tstep=1.0,
                   individual=False,
                   top=None,
                   dtype='dataset',
                   frame_indices=None):
    '''calculate diffusion

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
    act = CpptrajActions.Action_Diffusion()
    dslist = CpptrajDatasetList()

    _tsep = 'time ' + str(tstep)
    _individual = 'individual' if individual else ''

    # add 'df' as label
    label = 'df'
    command = ' '.join((mask, label, _tsep, _individual))

    # normally we just need
    # act(command, traj, top=_top, dslist=dslist)
    # but cpptraj need correct frame idx

    act.read_input(command, top=top, dslist=dslist)
    act.process(top)
    for idx, frame in enumerate(traj):
        # do not need mass
        act.do_action(frame, idx=idx)
    act.post_process()

    # make the label nicer
    for d in dslist:
        d.key = d.key.replace('[', '').replace(']', '').replace(label, '')

    return _get_data_from_dtype(dslist, dtype=dtype)


@_register_pmap
@_register_openmp
def calc_watershell(traj=None,
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

    Notes
    -----
    This method is not validated with cpptraj's output yet

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
    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)

    dslist = CpptrajDatasetList()

    act = CpptrajActions.Action_Watershell()

    if _solutemask in [None, '']:
        raise ValueError('must provide solute mask')

    _solventmask = solvent_mask if solvent_mask is not None else ''
    _noimage = 'noimage' if not image else ''

    _lower = 'lower ' + str(lower)
    _upper = 'upper ' + str(upper)
    command = ' '.join((_solutemask, _lower, _upper, _noimage, _solventmask))

    act(command, fi, top=_top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)


@_super_dispatch()
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
    >>> from pytraj.common_actions import calc_matrix
    >>> traj = pt.datafiles.load_trpcage()
    >>> mat = calc_matrix(traj, 'covar @CA')
    >>> # this is equal to
    >>> mat2 = pt.matrix.covar(traj, '@CA')
    >>> import numpy as np
    >>> np.testing.assert_equal(mat, mat2)
    '''
    act = CpptrajActions.Action_Matrix()

    dslist = CpptrajDatasetList()
    act(mask, traj, top=top, dslist=dslist)
    act.post_process()
    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
@_super_dispatch()
def calc_radgyr(traj=None,
                mask="",
                top=None,
                nomax=True,
                frame_indices=None,
                dtype='ndarray',
                *args,
                **kwd):
    '''calc radgyr

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.radgyr(traj, '@CA')
    >>> data = pt.radgyr(traj, '!:WAT', nomax=False)
    >>> data = pt.radgyr(traj, '@CA', frame_indices=[2, 4, 6])
    '''
    _nomax = 'nomax' if nomax else ""
    command = " ".join((mask, _nomax))

    act = CpptrajActions.Action_Radgyr()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
@_super_dispatch()
def calc_molsurf(traj=None,
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
    _offset = 'offset ' + str(offset) if offset != 0. else ''
    command = ' '.join((mask, _probe, _offset))

    act = CpptrajActions.Action_Molsurf()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
@_super_dispatch(has_ref=True)
def calc_rotation_matrix(traj=None,
                         ref=0,
                         mask="",
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
    dslist = CpptrajDatasetList()
    _mass = 'mass' if mass else ''

    command = ' '.join(('tmp', mask, 'savematrices', _mass))

    act = CpptrajActions.Action_Rmsd()
    act(command, [ref, traj], top=top, dslist=dslist)
    mat = dslist[-1].values
    # exclude data for reference
    if with_rmsd:
        return mat[1:], np.array(dslist[0].values[1:])
    else:
        return mat[1:]


@_super_dispatch()
def calc_volume(traj=None,
                mask="",
                top=None,
                dtype='ndarray',
                frame_indices=None):
    '''calculate volume

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> vol = pt.calc_volume(traj, '@CA')
    '''
    command = mask

    act = CpptrajActions.Action_Volume()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype)


@_super_dispatch()
def calc_multivector(traj,
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
    >>> vecs = pt.calc_multivector(traj, resrange='1-5', names=('C', 'N'))
    >>> vecs = pt.calc_multivector(traj, resrange='1-5', names='C N')
    '''
    act = CpptrajActions.Action_MultiVector()

    _resrange = 'resrange ' + resrange
    if 'name1' in names or 'name2' in names:
        # cpptraj style
        _names = names
    else:
        if isinstance(names, string_types):
            name1, name2 = names.split()
        elif isinstance(names, (list, tuple)):
            name1, name2 = names
        else:
            name1, name2 = '', ''
        _names = ' '.join(('name1', name1, 'name2', name2))
    command = ' '.join((_resrange, _names))

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
@_super_dispatch()
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

    _grid_spacing = ' '.join([str(x) for x in grid_spacing])
    _radscale = 'radscale ' + str(radscale)
    _buffer = 'buffer ' + str(buffer)
    _peakcut = 'peakcut ' + str(peakcut)
    _centermask = 'centermask ' + centermask

    if isinstance(size, tuple):
        assert len(size) == 3, 'lenghth of size must be 3'
    elif size is not None:
        raise ValueError('size must be None or a tuple. Please check method doc')

    _size = '' if size is None else 'size ' + ','.join([str(x) for x in size])

    if _size:
        # ignore buffer
        _buffer = ''
        # add center
        if center is not None:
            _center = 'center ' + ','.join([str(x) for x in center])
        else:
            _center = ''
    else:
        _center = ''

    command = ' '.join((dummy_filename, _grid_spacing, _center, _size, mask, _radscale, _buffer,
                        _centermask, _peakcut))

    act = CpptrajActions.Action_Volmap()

    dslist = CpptrajDatasetList()
    act(command, traj, top=top, dslist=dslist)
    act.post_process()
    return _get_data_from_dtype(dslist, dtype)


calc_volmap = volmap


@_register_openmp
def calc_rdf(traj=None,
             solvent_mask=':WAT@O',
             solute_mask='',
             maximum=10.,
             bin_spacing=0.5,
             image=True,
             density=0.033456,
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

    traj = _get_fiterator(traj, frame_indices)
    _top = _get_topology(traj, top)

    act = CpptrajActions.Action_Radial()

    if not isinstance(solvent_mask, string_types):
        solvent_mask = array_to_cpptraj_atommask(solvent_mask)

    if not isinstance(solute_mask, string_types) and solute_mask is not None:
        solute_mask = array_to_cpptraj_atommask(solute_mask)

    _spacing = str(bin_spacing)
    _maximum = str(maximum)
    _solventmask = solvent_mask
    _solutemask = solute_mask
    _noimage = 'noimage' if not image else ''
    _density = 'density ' + str(density) if density is not None else ''
    _center1 = 'center1' if center_solvent else ''
    _center2 = 'center2' if center_solute else ''
    _nointramol = 'nointramol' if not intramol else ''

    # order does matters
    # the order between _solventmask and _solutemask is swapped compared
    # to cpptraj's doc (to get correct result)
    command = ' '.join(
        ("pytraj_tmp_output.agr", _spacing, _maximum, _solutemask,
         _solventmask, _noimage, _density, _center1, _center2, _nointramol))

    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist)
    act.post_process()

    # make a copy sine dslist[-1].values return view of its data
    # dslist will be freed
    values = np.array(dslist[-1].values)
    # return (bin_centers, values)
    return (np.arange(bin_spacing / 2., maximum, bin_spacing), values)


@_super_dispatch()
def calc_pairdist(traj,
                  mask="*",
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
    dslist = CpptrajDatasetList()
    act = CpptrajActions.Action_PairDist()

    _command = 'mask ' + mask
    _delta = 'delta ' + str(delta)
    command = ' '.join((_command, _delta))

    act(command, traj, top=top, dslist=dslist)
    act.post_process()

    return _get_data_from_dtype(dslist, dtype=dtype)


pairdist = calc_pairdist


@_super_dispatch()
def calc_jcoupling(traj=None,
                   mask="",
                   top=None,
                   kfile=None,
                   dtype='dataset',
                   frame_indices=None):
    """
    Parameters
    ----------
    traj : any things that make `frame_iter_master` returning Frame object
    command : str, default ""
        cpptraj's command/mask
    kfile : str, default None, optional
        Dir for Karplus file. If "None", use $AMBERHOME dir
    dtype : str, {'dataset', ...}, default 'dataset'
    *args, **kwd: optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.calc_jcoupling(traj, ':1-12', kfile='data/Karplus.txt')
    """
    command = mask

    act = CpptrajActions.Action_Jcoupling()
    dslist = CpptrajDatasetList()

    if kfile is not None:
        command += " kfile %s" % kfile
    act(command, traj, dslist=dslist, top=top)
    return _get_data_from_dtype(dslist, dtype)


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

    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)

    CpptrajActions.Action_Translate()(command, fi, top=_top)


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
    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)
    CpptrajActions.Action_Scale()(command, fi, top=_top)


scale = do_scaling


def do_rotation(traj=None, command="", frame_indices=None, top=None):
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
    _top = _get_topology(traj, top)
    _assert_mutable(traj)
    fi = _get_fiterator(traj, frame_indices)
    CpptrajActions.Action_Rotate()(command, fi, top=_top)


rotate = do_rotation


@_super_dispatch()
def do_autoimage(traj, command="", frame_indices=None, top=None):
    '''perform autoimage and return the coordinate-updated traj

    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()[:]
    >>> traj = pt.autoimage(traj)
    '''
    _assert_mutable(traj)
    CpptrajActions.Action_AutoImage()(command, traj, top=top)
    return traj


autoimage = do_autoimage


@_register_pmap
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
    _top = _get_topology(traj, top)
    try:
        fi = traj.iterframe(autoimage=autoimage,
                            rmsfit=rmsfit,
                            frame_indices=frame_indices)
    except AttributeError:
        fi = _get_fiterator(traj, frame_indices)

    dslist = CpptrajDatasetList()
    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)
    else:
        mask = mask

    # add "crdset s1" to trick cpptraj dump coords to DatSetList
    command = mask + " crdset s1"

    act = CpptrajActions.Action_Average()
    act(command, fi, _top, dslist=dslist)

    # need to call this method so cpptraj will write
    act.post_process()

    frame = dslist[0].get_frame()
    if dtype.lower() == 'frame':
        return frame
    elif dtype.lower() in ['traj', 'trajectory']:
        new_top = _top if mask is '' else _top[mask]
        return Trajectory(xyz=frame.xyz.reshape(1, frame.n_atoms, 3).copy(),
                          top=new_top)
    else:
        raise ValueError('dtype must be frame or trajectory')


get_average_frame = mean_structure


@_register_pmap
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


@_super_dispatch()
def randomize_ions(traj=None, mask="", top=None, frame_indices=None):
    """randomize_ions for given Frame with Topology

    Parameters
    ----------
    traj : Trajectory-like or a Frame
        ``traj`` must be mutable
    mask : str
        cpptraj command
    frame_indices : {None, array-like}, optional
    top : Topology, optional (only needed if ``traj`` does not have it)
    """
    command = mask
    _assert_mutable(traj)

    act = CpptrajActions.Action_RandomizeIons()
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
    from pytraj.analyses.CpptrajAnalyses import Analysis_Clustering

    dslist = CpptrajDatasetList()
    dslist.add_set('double', '__array_like')
    dslist[0].resize(len(array_like))
    dslist[0].values[:] = array_like
    act = Analysis_Clustering()
    command = 'data __array_like ' + command
    act(command, dslist=dslist)

    return np.array(dslist[-1])


@_register_pmap
@_super_dispatch()
def calc_multidihedral(traj=None,
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
    *arg and **kwd: additional arguments (for advanced users)


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
        if is_int(resrange):
            resrange = [resrange, ]
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

    _command = " ".join((d_types, _resrange, dh_types, _range360))

    dslist = CpptrajDatasetList()
    act = CpptrajActions.Action_MultiDihedral()
    act(_command, traj, top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)


@_super_dispatch()
def calc_atomicfluct(traj=None,
                     mask="",
                     top=None,
                     dtype='ndarray',
                     frame_indices=None):
    '''

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.calc_atomicfluct(traj, '@CA')
    >>> data[:3]
    array([[  5.        ,   0.61822273],
           [ 16.        ,   0.5627449 ],
           [ 40.        ,   0.53717119]])
    '''
    dslist = CpptrajDatasetList()
    act = CpptrajActions.Action_AtomicFluct()
    act(mask, traj, top=top, dslist=dslist)
    act.post_process()
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_bfactors(traj=None,
                  mask="",
                  byres=True,
                  top=None,
                  dtype='ndarray',
                  frame_indices=None):
    # Not: do not use _super_dispatch here since we used in calc_atomicfluct
    """calculate pseudo bfactor

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
    # do not use _super_dispatch again
    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)
    _command = " ".join((mask, byres_text, "bfactor"))
    return calc_atomicfluct(traj=traj,
                            mask=_command,
                            top=top,
                            dtype=dtype,
                            frame_indices=frame_indices)


@_register_pmap
def calc_vector(traj=None,
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
    from pytraj.core.action_list import ActionList

    dslist = CpptrajDatasetList()
    _top = _get_topology(traj, top)
    list_of_commands = _get_list_of_commands(command)
    fi = _get_fiterator(traj, frame_indices)
    actlist = ActionList()

    for command in list_of_commands:
        act = CpptrajActions.Action_Vector()
        actlist.add_action(act, command, _top, dslist=dslist)
    actlist.do_actions(fi)

    return _get_data_from_dtype(dslist, dtype=dtype)


@_super_dispatch()
def _calc_vector_center(traj=None,
                        mask="",
                        top=None,
                        mass=False,
                        dtype='ndarray',
                        frame_indices=None):

    dslist = CpptrajDatasetList()
    dslist.set_own_memory(False)  # need this to avoid segmentation fault
    act = CpptrajActions.Action_Vector()
    command = "center " + mask

    if mass:
        command += " mass"

    act(command, traj, top=top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)


@_register_pmap
def calc_center_of_mass(traj=None,
                        mask='',
                        top=None,
                        dtype='ndarray',
                        frame_indices=None):
    # note: do not use _super_dispatch for this method since
    # we already use for _calc_vector_center
    return _calc_vector_center(traj=traj,
                               mask=mask,
                               top=top,
                               mass=True,
                               dtype=dtype,
                               frame_indices=frame_indices)


calc_COM = calc_center_of_mass


@_register_pmap
@_super_dispatch()
def calc_center_of_geometry(traj=None,
                            mask="",
                            top=None,
                            dtype='ndarray',
                            frame_indices=None):

    atom_mask_obj = top(mask)
    dslist = CpptrajDatasetList()
    dslist.add_set("vector")

    for frame in iterframe_master(traj):
        dslist[0].append(frame.center_of_geometry(atom_mask_obj))
    return _get_data_from_dtype(dslist, dtype=dtype)


calc_COG = calc_center_of_geometry


# do not use _super_dispatch here since we did in inside this method
# to avoid complicated code checking.

@_register_openmp
def calc_pairwise_rmsd(traj=None,
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
    *args, **kwd: optional (for advanced user)

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

    Notes
    -----
    Install ``libcpptraj`` with ``openmp`` to get benefit from parallel
    """
    # we copy Frame coordinates to DatasetCoordsCRD first

    if not isinstance(mask, string_types):
        mask = array_to_cpptraj_atommask(mask)

    act = CpptrajAnalyses.Analysis_Rms2d()

    crdname = 'default_coords'
    dslist, _top, command = _get_fi_with_dslist(traj, mask, frame_indices, top, crdname=crdname)

    command = ' '.join((command, metric, "crdset {} rmsout mycrazyoutput".format(crdname)))

    act(command, _top, dslist=dslist)
    # remove dataset coords to free memory
    dslist.remove_set(dslist[0])

    if dtype == 'ndarray':
        return _get_matrix_from_dataset(dslist[0], mat_type)
    else:
        return _get_data_from_dtype(dslist, dtype)



@_register_pmap
@_super_dispatch()
def calc_temperatures(traj=None,
                      mask="",
                      frame_indices=None,
                      top=None,
                      dtype='ndarray'):
    """return 1D python array of temperatures (from velocity) in traj
    if `frame` keyword is specified cpptraj/pytraj will take existing T

    Default = array of 0.0
    """
    dslist = CpptrajDatasetList()

    act = CpptrajActions.Action_Temperature()
    act(mask, traj, dslist=dslist, top=top)

    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
def rmsd_perres(traj=None,
                ref=0,
                mask="",
                mass=False,
                resrange=None,
                perres_mask=None,
                perres_center=False,
                perres_invert=False,
                frame_indices=None,
                top=None,
                dtype='dataset'):
    """superpose ``traj`` to ``ref`` with `mask`, then calculate nofit rms for residues
    in `resrange` with given `perresmask`

    Returns
    -------
    out : pytraj.DatasetList, shape=(1+n_residues, n_frames)
        out[0]: regular rmsd
        out[1:]: perres rmsd for all given residues
        `out.values` will return corresponding numpy array
    """
    _range = 'range %s ' % resrange
    _perresmask = 'perresmask ' + perres_mask if perres_mask is not None else ''
    _perrestcenter = 'perrescenter' if perres_center else ''
    _perrestinvert = 'perresinvert' if perres_invert else ''

    cm = " ".join((mask, 'perres', _range, _perresmask, _perrestcenter,
                   _perrestinvert))
    return calc_rmsd(traj=traj,
                     ref=ref,
                     mask=cm,
                     nofit=False,
                     mass=mass,
                     frame_indices=frame_indices,
                     top=top,
                     dtype=dtype)


@_register_pmap
def calc_rmsd_nofit(traj=None,
                    ref=0,
                    mask="",
                    mass=False,
                    frame_indices=None,
                    top=None,
                    dtype='ndarray'):
    '''
    See also
    --------
    calc_rmsd
    '''
    return calc_rmsd(traj=traj,
                     ref=ref,
                     mask=mask,
                     mass=mass,
                     nofit=True,
                     frame_indices=frame_indices,
                     top=top,
                     dtype=dtype)


@_register_pmap
def calc_rmsd(traj=None,
              ref=0,
              mask="",
              nofit=False,
              mass=False,
              frame_indices=None,
              top=None,
              dtype='ndarray'):
    """calculate rmsd

    Parameters
    ----------
    traj : Trajectory-like
    ref : {Frame, int}, default=0 (first frame)
        Reference frame or index.
    mask : str or 1D array-like of string or 1D or 2D array-like
        Atom mask/indices
    nofit : bool, default False
        if False, perform fitting (rotation and translation).
        if ``traj`` is mutable, its coordinates will be updated
        if True, not fitting.
    top : {Topology, str}, default None, optional
    dtype : return data type, default='ndarray'


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

    Notes
    -----
    if ``traj`` is mutable, its coordinates will be updated

    """

    _nofit = ' nofit ' if nofit else ''
    _mass = ' mass ' if mass else ''
    opt = _nofit + _mass

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
        elif 'int' in dname or 'object' in dname:
            if cmd.ndim == 1 and 'object' not in dname:
                command = [array_to_cpptraj_atommask(mask), ]
            elif cmd.ndim == 2 or 'object' in dname:
                command = [array_to_cpptraj_atommask(x) for x in mask]
            else:
                raise ValueError("only support array with ndim=1,2")
        else:
            raise ValueError("not supported")

    _top = _get_topology(traj, top)

    ref = _get_reference_from_traj(traj, ref)
    fi = _get_fiterator(traj, frame_indices)

    alist = ActionList()
    dslist = CpptrajDatasetList()

    for cm in command:
        _cm = cm + opt
        if 'savematrices' in _cm:
            if dtype not in ['dataset', 'cpptraj_dataset']:
                raise ValueError('if savematrices, dtype must be "dataset"')
            _cm = 'RMDSset ' + _cm
        alist.add_action(CpptrajActions.Action_Rmsd(),
                         _cm,
                         top=_top,
                         dslist=dslist)

    alist.do_actions(ref)
    alist.do_actions(fi)

    dnew = DatasetList(dslist)
    for d in dnew:
        d.values = d.values[1:]
    return _get_data_from_dtype(dnew, dtype=dtype)

# alias for `calc_rmsd`
rmsd = calc_rmsd


@_register_pmap
@_super_dispatch(has_ref=True)
def calc_distance_rmsd(traj=None,
                       ref=0,
                       mask='',
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
    dslist = CpptrajDatasetList()
    command = mask

    act = CpptrajActions.Action_DistRmsd()
    act(command, [ref, traj], top=top, dslist=dslist)

    # exclude ref value
    for d in dslist:
        d.data = d.data[1:]
    return _get_data_from_dtype(dslist, dtype=dtype)

# alias
distance_rmsd = calc_distance_rmsd


@_super_dispatch()
def align_principal_axis(traj=None, mask="*", top=None, frame_indices=None):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (Trajectory, Frame)
    """
    _assert_mutable(traj)

    act = CpptrajActions.Action_Principal()
    command = mask + " dorotation"
    act(command, traj, top=top)


def principal_axes(traj=None, mask='*', dorotation=False, mass=True, top=None):
    """calculate eigenvalues and eigenvectors

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, default '*' (all atoms)
    mass: bool, defaul True
    if `dorotation`, the system will be aligned along principal axes (apply for mutable system)
    top : Topology, optional

    Returns
    -------
    out_0 : numpy array, shape=(n_frames, 3, 3), corresponding eigenvectors
    out_1: numpy array with shape=(n_frames, 3), eigenvalues
    """
    act = CpptrajActions.Action_Principal()
    command = mask

    _dorotation = 'dorotation' if dorotation else ''
    _mass = 'mass' if mass else ''

    if 'name' not in command:
        command += ' name pout'

    command = ' '.join((command, _dorotation, _mass))

    _top = _get_topology(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, _top, dslist=dslist)
    return (dslist[0].values, dslist[1].values)


def _closest_iter(act, traj):
    '''

    Parameters
    ----------
    act : Action object
    traj : Trajectory-like
    '''

    for frame in iterframe_master(traj):
        new_frame = act.do_action(frame, get_new_frame=True)
        yield new_frame


@_register_openmp
@_super_dispatch()
def closest(traj=None,
            mask='*',
            solvent_mask=None,
            n_solvents=10,
            frame_indices=None,
            top=None):
    """return either a new Trajectory or a frame iterator. Keep only ``n_solvents`` closest to mask

    Parameters
    ----------
    traj : Trajectory-like | list of Trajectory-like/frames | frame iterator | chunk iterator
    mask: str, default '*' (all solute atoms)
    top : Topology-like object, default=None, optional

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
    dslist = CpptrajDatasetList()

    command = str(n_solvents) + ' ' + mask

    act = CpptrajActions.Action_Closest()

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

    act.read_input(command, top, dslist=dslist)
    new_top = act.process(top, get_new_top=True)

    fiter = _closest_iter(act, traj)

    return (fiter, new_top.copy())


@_register_pmap
@_super_dispatch(has_ref=True)
def native_contacts(traj=None,
                    ref=0,
                    mask="",
                    mask2="",
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
    act = CpptrajActions.Action_NativeContacts()
    dslist = CpptrajDatasetList()

    if not isinstance(mask2, string_types):
        # [1, 3, 5] to "@1,3,5
        mask2 = array_to_cpptraj_atommask(mask2)
    command = ' '.join((mask, mask2))

    _distance = str(distance)
    _noimage = "noimage" if not image else ""
    _includesolvent = "includesolvent" if include_solvent else ""
    _byres = "byresidue" if byres else ""

    _command = " ".join(('ref myframe', command, _distance, _noimage,
                         _includesolvent, _byres))
    dslist.add_set('ref_frame', 'myframe')
    dslist[0].add_frame(ref)
    dslist[0].top = top
    act(_command, traj, top=top, dslist=dslist)
    dslist._pop(0)

    return _get_data_from_dtype(dslist, dtype=dtype)


@_super_dispatch()
def calc_grid(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """
    """
    # TODO: doc, rename method, move to seperate module?
    act = CpptrajActions.Action_Grid()
    dslist = CpptrajDatasetList()

    # cpptraj require output
    command = "tmp_pytraj_grid_output.txt " + command
    with goto_temp_folder():
        act(command, traj, dslist=dslist, top=top, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)


def check_structure(traj=None, command="", top=None, *args, **kwd):
    """check if the structure is ok or not

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_rna()
    >>> check_structure(traj[0], top=traj.top)
    """
    act = CpptrajActions.Action_CheckStructure()

    # cpptraj require output
    _top = _get_topology(traj, top)
    act(command, traj, top=_top, *args, **kwd)


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
    act = CpptrajAnalyses.Analysis_Timecorr()

    cdslist = CpptrajDatasetList()

    cdslist.add_set("vector", "_vec0")
    cdslist.add_set("vector", "_vec1")
    cdslist[0].data = np.asarray(vec0).astype('f8')
    cdslist[1].data = np.asarray(vec1).astype('f8')

    _order = "order " + str(order)
    _tstep = "tstep " + str(tstep)
    _tcorr = "tcorr " + str(tcorr)
    _norm = "norm" if norm else ""
    command = " ".join(('vec1 _vec0 vec2 _vec1', _order, _tstep, _tcorr, _norm
                        ))
    act(command, dslist=cdslist)
    return _get_data_from_dtype(cdslist[2:], dtype=dtype)


def crank(data0, data1, mode='distance', dtype='ndarray'):
    """
    Parameters
    ----------
    data0 : array-like
    data1 : array-like
    mode : str, {'distance', 'angle'}
    dtype : str

    Notes
    -----
    Same as `crank` in cpptraj
    """
    from pytraj.analyses.CpptrajAnalyses import Analysis_CrankShaft

    cdslist = CpptrajDatasetList()
    cdslist.add_set("double", "d0")
    cdslist.add_set("double", "d1")

    cdslist[0].data = np.asarray(data0)
    cdslist[1].data = np.asarray(data1)

    act = Analysis_CrankShaft()
    command = ' '.join((mode, 'd0', 'd1'))
    act(command, dslist=cdslist)
    return _get_data_from_dtype(cdslist[2:], dtype=dtype)


def cross_correlation_function(data0, data1, dtype='ndarray'):
    """
    Notes
    -----
    Same as `corr` in cpptraj
    """

    cdslist = CpptrajDatasetList()
    cdslist.add_set("double", "d0")
    cdslist.add_set("double", "d1")

    cdslist[0].data = np.asarray(data0)
    cdslist[1].data = np.asarray(data1)

    act = CpptrajAnalyses.Analysis_Corr()
    act("d0 d1 out _tmp.out", dslist=cdslist)
    return _get_data_from_dtype(cdslist[2:], dtype=dtype)


def auto_correlation_function(data, dtype='ndarray', covar=True):
    """
    Notes
    -----
    Same as `autocorr` in cpptraj
    """

    _nocovar = " " if covar else " nocovar"

    cdslist = CpptrajDatasetList()
    cdslist.add_set("double", "d0")

    cdslist[0].data = np.asarray(data)

    act = CpptrajAnalyses.Analysis_AutoCorr()
    command = "d0 out _tmp.out" + _nocovar
    act(command, dslist=cdslist)
    return _get_data_from_dtype(cdslist[1:], dtype=dtype)


def lifetime(data, cut=0.5, rawcurve=False, more_options='', dtype='ndarray'):
    """lifetime (adapted lightly from cpptraj doc)

    Parameters
    ----------
    data : 1D-array or 2D array-like
    cut : cutoff to use when determining if data is 'present', default 0.5
    more_options : str, more cpptraj's options. Check cpptraj's manual.
    """
    data = np.asarray(data)
    if data.ndim == 1:
        data_ = [data, ]
    else:
        data_ = data

    _outname = 'name lifetime_'
    _cut = 'cut ' + str(cut)
    _rawcurve = 'rawcurve' if rawcurve else ''
    # do not sorting dataset's names. We can accessing by indexing them.
    _nosort = 'nosort'

    namelist = []
    cdslist = CpptrajDatasetList()
    for idx, arr in enumerate(data_):
        # create datasetname so we can reference them
        name = 'data_' + str(idx)
        if 'int' in arr.dtype.name:
            cdslist.add_set("integer", name)
        else:
            cdslist.add_set("double", name)
        cdslist[-1].data = np.asarray(arr)
        namelist.append(name)

    act = CpptrajAnalyses.Analysis_Lifetime()
    _cm = ' '.join(namelist)
    command = " ".join((_cm, _outname, _cut, _rawcurve, _nosort, more_options))
    act(command, dslist=cdslist)

    for name in namelist:
        cdslist.remove_set(cdslist[name])
    return _get_data_from_dtype(cdslist, dtype=dtype)


@_super_dispatch()
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
    dslist = DatasetList()

    for idx, frame in enumerate(iterframe_master(traj)):
        top.set_distance_mask_reference(frame)
        dslist.append({str(idx): np.asarray(top.select(mask))})
    return _get_data_from_dtype(dslist, dtype)


@_register_pmap
def pucker(traj=None,
           pucker_mask=("C1'", "C2'", "C3'", "C4'", "O4'"),
           resrange=None,
           top=None,
           dtype='dataset',
           range360=False,
           method='altona',
           use_com=True,
           amplitude=True,
           offset=None,
           *args,
           **kwd):
    """Note: not validate yet

    """
    from pytraj.compat import range

    _top = _get_topology(traj, top)
    if not resrange:
        resrange = range(_top.n_residues)

    _range360 = "range360" if range360 else ""
    geom = "geom" if not use_com else ""
    amp = "amplitude" if amplitude else ""
    _offset = "offset " + str(offset) if offset else ""

    cdslist = CpptrajDatasetList()
    for res in resrange:
        act = CpptrajActions.Action_Pucker()
        command = " ".join((":" + str(res + 1) + '@' + x for x in pucker_mask))
        name = "pucker_res" + str(res + 1)
        command = " ".join((name, command, _range360, method, geom, amp,
                            _offset))
        act(command, traj, top=_top, dslist=cdslist, *args, **kwd)
    return _get_data_from_dtype(cdslist, dtype)


@_super_dispatch()
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
    _center = '' if center == 'box' else center
    _mass = 'mass' if mass else ''
    command = ' '.join((mask, _center, _mass))
    _assert_mutable(traj)

    act = CpptrajActions.Action_Center()
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
    _top = _get_topology(traj, top)

    if "custom:" in mask:
        command = mask
    else:
        command = "custom:" + mask

    act = CpptrajActions.Action_MakeStructure()

    act(command, traj, top=_top)
    return traj


@_register_openmp
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
    _top = _get_topology(traj, top)
    if isinstance(direction, string_types):
        _direction = direction
    elif isinstance(direction, (list, tuple)):
        # example: direction = ('001, '0-10')
        _direction = 'dir ' + ' dir '.join(direction)
    else:
        raise ValueError(
            'only support ``direction`` as a string or list/tuple of strings')
    command = ' '.join(('name tmp_cell', _direction, mask))

    act = CpptrajActions.Action_ReplicateCell()
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist)
    traj = Trajectory(xyz=dslist[0].xyz, top=dslist[0].top)

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
    _top = _get_topology(traj, top)

    if not isinstance(resid, string_types):
        resid = str(resid + 1)
    deg = str(deg)

    command = ':'.join((dihedral_type, resid, dihedral_type, deg))
    make_structure(traj, command, top=_top)
    return traj


def make_structure(traj=None, mask="", top=None):
    _assert_mutable(traj)
    _top = _get_topology(traj, top)

    command = mask
    act = CpptrajActions.Action_MakeStructure()
    act(command, traj, top=_top)


@_super_dispatch()
def _projection(traj,
                mask,
                eigenvalues,
                eigenvectors,
                scalar_type,
                average_coords=None,
                frame_indices=None,
                dtype='ndarray',
                top=None):

    act = CpptrajActions.Action_Projection()
    dslist = CpptrajDatasetList()

    mode_name = 'my_modes'
    dslist.add_set('modes', mode_name)
    is_reduced = False
    dataset_mode = dslist[-1]
    n_vectors = len(eigenvalues)
    dataset_mode._set_modes(is_reduced, n_vectors, eigenvectors.shape[1],
                          eigenvalues, eigenvectors.flatten())
    dataset_mode.scalar_type = scalar_type

    if average_coords is not None:
        dataset_mode._allocate_avgcoords(3*average_coords.shape[0])
        dataset_mode._set_avg_frame(average_coords.flatten())

    _mask = mask
    _evecs = 'evecs {}'.format(mode_name)
    _beg_end = 'beg 1 end {}'.format(n_vectors)
    command = ' '.join((_evecs, _mask, _beg_end))
    act(command, traj, top=top, dslist=dslist)
    dslist._pop(0)

    return _get_data_from_dtype(dslist, dtype=dtype)


@_super_dispatch(has_ref=True)
def superpose(traj, ref=0, mask='', frame_indices=None, top=None): 
    act = CpptrajActions.Action_Rmsd()
    act(mask, traj, top=top)
    return traj


def pca(traj,
        mask,
        n_vecs=2,
        dtype='ndarray',
        top=None):
    '''perform PCA analysis by following below steps:

    - perform rmsfit to first frame with ``mask``
    - compute average frame with ``mask``
    - rmsfit to average frame with ``mask``
    - compute covariance matrix
    - diagonalize the matrix to get eigenvectors and eigenvalues
    - perform projection of each frame with mask to each eigenvector 

    Parameters
    ----------
    traj : Trajectory
        traj must be ``pytraj.Trajectory``, which can be created by ``pytraj.load`` method.
    mask : str
        atom mask
    n_vecs : int, default 2
        number of eigenvectors. If user want to compute projection for all eigenvectors, 
        should specify n_vecs=-1 (or a negative number)
    dtype : return datatype
    top : Topology, optional

    Returns
    -------
    projection_data: ndarray, shape=(n_vecs, n_frames)
    (eigenvalues, eigenvectors)

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
    '''
    # TODO: move to another file
    # NOTE: do not need to use _super_dispatch here since we already use in _projection
    from pytraj import matrix

    traj.superpose(ref=0, mask=mask)
    avg = mean_structure(traj)
    traj.superpose(ref=avg, mask=mask)
    avg2 = mean_structure(traj, mask=mask)

    mat = matrix.covar(traj, mask)
    if n_vecs < 0:
        n_vecs = mat.shape[0]
    else:
        n_vecs = n_vecs

    eigenvalues, eigenvectors = matrix.diagonalize(mat, n_vecs=n_vecs, dtype='tuple')
    projection_data = _projection(traj, mask=mask, average_coords=avg2.xyz,
                                  eigenvalues=eigenvalues, 
                                  eigenvectors=eigenvectors,
                                  scalar_type='covar', dtype=dtype)
    return projection_data, (eigenvalues, eigenvectors)

calc_pca = pca


@_super_dispatch()
@_register_openmp
def calc_atomiccorr(traj,
                    mask='',
                    cut=None,
                    min_spacing=None,
                    byres=True,
                    frame_indices=None,
                    dtype='ndarray',
                    top=None):
    '''Calculate average correlations between the motion of atoms in mask.

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

    act = CpptrajActions.Action_AtomicCorr()
    dslist = CpptrajDatasetList()

    with goto_temp_folder():
        act(command, traj, top=top, dslist=dslist)
        # need to post_process for this Action
        act.post_process()

    return _get_data_from_dtype(dslist, dtype=dtype)


def _rotdif(arr,
            nvecs=1000,
            rvecin=None,
            rseed=80531,
            order=2,
            ncorr=-1,
            tol=1E-6,
            d0=0.03,
            nmesh=2,
            dt=0.002,
            ti=0.0,
            tf=-1,
            itmax=-1.,
            dtype='ndarray'):
    '''

    Parameters
    ----------
    arr : array-like, shape (n_frames, 3, 3) or (n_frames, 9)
    '''
    _nvecs = 'nvecs ' + str(nvecs)
    _rvecin = 'rvecin ' + rvecin if rvecin is not None else ''
    _rseed = 'rseed ' + str(rseed)
    _order = 'order ' + str(order)
    _ncorr = 'ncorr ' + str(ncorr) if ncorr > 0 else ''
    _tol = 'tol ' + str(tol)
    _d0 = 'd0 ' + str(d0)
    _nmesh = 'nmesh ' + str(nmesh)
    _dt = 'dt ' + str(dt)
    _ti = 'ti ' + str(ti)
    _tf = 'tf ' + str(tf)
    _itmax = 'itmax ' + str(itmax)
    _rmatrix = 'rmatrix mymat'

    act = CpptrajAnalyses.Analysis_Rotdif()
    dslist = CpptrajDatasetList()
    dslist.add_set('mat3x3', 'mymat')

    arr = np.asarray(arr, dtype='f8')
    msg = 'array must have shape=(n_frames, 9) or (n_frames, 3, 3)'
    shape = arr.shape
    if arr.ndim == 2:
        assert shape[1] == 9, msg
    elif arr.ndim == 3:
        assert shape[1:3] == (3, 3), msg
        # need to reshape to (n_frames, 9)
        arr = arr.reshape(shape[0], shape[1] * shape[2])
    else:
        raise ValueError(msg)

    dslist[0]._append_from_array(arr)

    command = ' '.join((_nvecs, _rvecin, _rseed, _order, _ncorr, _tol, _d0,
                        _nmesh, _dt, _ti, _tf, _itmax, _rmatrix))

    act(command, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)


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
    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)
    act = CpptrajActions.Action_Bounds()
    dslist = CpptrajDatasetList()
    dx, dy, dz = grid_spacing
    _dx = 'dx ' + str(dx) if dx > 0. else ''
    _dy = 'dy ' + str(dy) if dy > 0. else ''
    _dz = 'dz ' + str(dz) if dz > 0. else ''
    _offset = 'offset ' + str(offset)
    command = ' '.join((mask, 'out tmp_bounds.dat', _dx, _dy, _dz,
                        'name grid_', _offset))

    with goto_temp_folder():
        act(command, fi, top=_top, dslist=dslist)
    act.post_process()

    return _get_data_from_dtype(dslist, dtype=dtype)


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
    '''calculate lowest curve for data

    Paramters
    ---------
    data : 2D array-like
    points : number of lowest points in each bin, default 10
    step : step size, default 0.2

    Return
    ------
    2d array
    '''
    _points = 'points ' + str(points)
    _step = 'step ' + str(step)
    label = 'mydata'
    command = ' '.join((label, _points, _step))

    data = np.asarray(data)

    act = CpptrajAnalyses.Analysis_LowestCurve()
    dslist = CpptrajDatasetList()

    dslist.add_new('xymesh', label)
    dslist[0]._append_from_array(data.T)

    act(command, dslist=dslist)
    return np.array([dslist[-1]._xcrd(), np.array(dslist[-1].values)])
