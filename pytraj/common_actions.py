"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd
>>> from pytraj.common_actions import translate
# TODO : use __all__
"""
from __future__ import absolute_import
import os
from array import array
from functools import partial

from pytraj.action_dict import ActionDict
adict = ActionDict()

from pytraj.analysis_dict import AnalysisDict
analdict = AnalysisDict()

from ._get_common_objects import _get_top, _get_data_from_dtype, _get_list_of_commands
from ._get_common_objects import _get_matrix_from_dataset
from ._get_common_objects import _get_reference_from_traj
from ._common_actions import calculate
from .utils import _import_numpy, is_array, ensure_not_none_or_string
from .utils import is_int
from .utils.context import goto_temp_folder
from .utils.convert import array_to_cpptraj_atommask as to_cpptraj_atommask
from .externals.six import string_types
from .Frame import Frame
#from .Trajectory import Trajectory
from .AtomMask import AtomMask
from .Topology import Topology
from .datasets.DataSetList import DataSetList as CpptrajDatasetList
from .core.DataFileList import DataFileList
from .math.DistRoutines import distance
from .externals.gdt.calc_score import calc_score
from .hbonds import search_hbonds, search_nointramol_hbonds
from .dssp_analysis import calc_dssp
from ._nastruct import nastruct
from ._shared_methods import _frame_iter_master
from .externals.get_pysander_energies import get_pysander_energies
from . import _long_manual
from .decorators import noparallel

list_of_cal = ['calc_distance',
               'calc_dihedral',
               'calc_radgyr',
               'calc_angle',
               'calc_molsurf',
               'calc_distrmsd',
               'calc_volume',
               'calc_protein_score',
               'calc_dssp',
               'calc_matrix',
               'calc_jcoupling',
               'calc_radial',
               'calc_watershell',
               'calc_vector',
               'calc_multivector',
               'calc_volmap',
               'calc_rdf',
               'calc_multidihedral',
               'calc_atomicfluct',
               'calc_COM',
               'calc_center_of_mass',
               'calc_center_of_geometry',
               'calc_pairwise_rmsd',
               'calc_density',
               'calc_grid',
               'calc_temperatures',
               'calc_linear_interaction_energy', ]

list_of_do = ['do_translation',
              'do_rotation',
              'do_autoimage',
              'do_clustering', ]

list_of_get = ['get_average_frame']

list_of_the_rest = ['search_hbonds', 'search_nointramol_hbonds',
                    'align_principal_axis', 'principal_axes', 'closest',
                    'native_contacts', 'nastruct']

__all__ = list_of_do + list_of_cal + list_of_get + list_of_the_rest

calc_protein_score = calc_score
calc_energies = get_pysander_energies
energy_decomposition = get_pysander_energies

action_type = calculate


def _noaction_with_TrajectoryIterator(trajiter):
    from pytraj import TrajectoryIterator
    if isinstance(trajiter, TrajectoryIterator):
        raise ValueError(
            "This analysis does not support immutable object. Use `pytraj.Trajectory`")


def calc_distance(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """calculate distance between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate distance for two atoms, using amber mask
    >>> pt.distance(traj, '@1 @3')

    >>> # calculate distance for two groups of atoms, using amber mask
    >>> pt.distance(traj, '@1,37,8 @2,4,6')

    >>> # calculate distance between two residues, using amber mask
    >>> pt.distance(traj, ':1 :10')

    >>> # calculate multiple distances between two residues, using amber mask
    >>> # distance between residue 1 and 10, distance between residue 3 and 20 
    >>> # (when using atom string mask, index starts from 1)
    >>> pt.distance(traj, [':1 :10', ':3 :20'])

    >>> # calculate distance for a series of atoms, using array for atom mask
    >>> # distance between atom 1 and 5, distance between atom 4 and 10 (index starts from 0)
    >>> pt.distance(traj, [[1, 5], [4, 10]])
    """
    import numpy as np
    ensure_not_none_or_string(traj)
    command = mask

    _top = _get_top(traj, top)

    cm_arr = np.asarray(command)

    if 'int' in cm_arr.dtype.name:
        from pytraj.datasetlist import from_dict

        int_2darr = cm_arr

        if int_2darr.ndim == 1:
            int_2darr = np.atleast_2d(cm_arr)

        if int_2darr.shape[1] != 2:
            raise ValueError("require int-array with shape=(n_atoms, 2)")

        if 'n_frames' not in kwd.keys():
            try:
                n_frames = traj.n_frames
            except:
                raise ValueError("require specifying n_frames")
        else:
            n_frames = kwd['n_frames']

        arr = np.empty([n_frames, len(int_2darr)])

        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_distance(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = from_dict({'distance': arr})
            return _get_data_from_dtype(py_dslist, dtype)

    elif isinstance(command, (list, tuple, string_types, np.ndarray)):
        # create a list
        if not isinstance(command, np.ndarray):
            list_of_commands = _get_list_of_commands(command)
        else:
            list_of_commands = command

        from pytraj.core.ActionList import ActionList
        from pytraj.actions.CpptrajActions import Action_Distance

        dslist = CpptrajDatasetList()
        actlist = ActionList()

        for cm in list_of_commands:
            actlist.add_action(
                Action_Distance(), cm, _top,
                dslist=dslist, *args, **kwd)
        actlist.do_actions(traj)
        return _get_data_from_dtype(dslist, dtype)

    else:
        raise ValueError(
            "command must be a string, a list/tuple of strings, or "
            "a numpy 2D array")


def calc_angle(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """calculate angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate angle for three atoms, using amber mask
    >>> pt.angle(traj, '@1 @3 @10')

    >>> # calculate angle for three groups of atoms, using amber mask
    >>> pt.angle(traj, '@1,37,8 @2,4,6 @5,20')

    >>> # calculate angle between three residues, using amber mask
    >>> pt.angle(traj, ':1 :10 :20')

    >>> # calculate multiple angles between three residues, using amber mask
    >>> # angle between residue 1, 10, 20, angle between residue 3, 20, 30
    >>> # (when using atom string mask, index starts from 1)
    >>> pt.angle(traj, [':1 :10 :20', ':3 :20 :30'])

    >>> # calculate angle for a series of atoms, using array for atom mask
    >>> # angle between atom 1, 5, 8, distance between atom 4, 10, 20 (index starts from 0)
    >>> pt.angle(traj, [[1, 5, 8], [4, 10, 20]])
    """
    import numpy as np
    from pytraj.datasetlist import from_dict
    command = mask

    ensure_not_none_or_string(traj)

    _, np = _import_numpy()
    _top = _get_top(traj, top)
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
            dset = calculate(
                "angle", traj, command,
                top=_top,
                quick_get=True, *args, **kwd)
            return _get_data_from_dtype(dset, dtype)

        elif isinstance(command, (list, tuple)):
            list_of_commands = command
            from pytraj.core.ActionList import ActionList
            from pytraj.actions.CpptrajActions import Action_Angle
            dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add_action(
                    Action_Angle(), cm, _top,
                    dslist=dslist, *args, **kwd)
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
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_angle(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            py_dslist = from_dict({'angle': arr})
            return _get_data_from_dtype(py_dslist, dtype)


def _dihedral_res(traj, mask=(), resid=0, dtype='ndarray', top=None):
    '''calculate dihedral within a single residue. For internal use only.

    Parameters
    ----------
    traj
    mask : tuple of strings
    resid
    dtype
    '''

    if is_int(resid):
        resid = str(resid + 1)
    else:
        resid = resid
    m = ' :%s@' % resid
    command = m + m.join(mask)
    return calc_dihedral(traj=traj, mask=command, top=top, dtype=dtype)


def calc_dihedral(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """calculate dihedral angle between two maskes

    Parameters
    ----------
    traj : Trajectory-like, list of Trajectory, list of Frames
    mask : str or array
    top : Topology, optional
    dtype : return type, defaul 'ndarray'

    Examples
    --------
    >>> import pytraj as pt
    >>> # calculate dihedral angle for four atoms, using amber mask
    >>> pt.dihedral(traj, '@1 @3 @10 @20')

    >>> # calculate dihedral angle for four groups of atoms, using amber mask
    >>> pt.dihedral(traj, '@1,37,8 @2,4,6 @5,20 @21,22')

    >>> # calculate dihedral angle for four residues, using amber mask
    >>> pt.dihedral(traj, ':1 :10 :20 :22')

    >>> # calculate multiple dihedral angles for four residues, using amber mask
    >>> # angle for residue 1, 10, 20, 30; angle between residue 3, 20, 30, 40
    >>> # (when using atom string mask, index starts from 1)
    >>> pt.dihedral(traj, [':1 :10 :20 :30', ':3 :20 :30 :40'])

    >>> # calculate dihedral angle for a series of atoms, using array for atom mask
    >>> # dihedral angle for atom 1, 5, 8, 10, dihedral for atom 4, 10, 20, 30 (index starts from 0)
    >>> pt.dihedral(traj, [[1, 5, 8, 10], [4, 10, 20, 30]])
    """
    import numpy as np
    ensure_not_none_or_string(traj)
    command = mask

    _, np = _import_numpy()
    _top = _get_top(traj, top)
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
            dset = calculate(
                "dihedral", traj, command,
                top=_top,
                quick_get=True, *args, **kwd)
            return _get_data_from_dtype(dset, dtype)

        elif isinstance(command, (list, tuple)):
            list_of_commands = command
            from pytraj.core.ActionList import ActionList
            from pytraj.actions.CpptrajActions import Action_Dihedral
            dslist = CpptrajDatasetList()
            actlist = ActionList()

            for cm in list_of_commands:
                actlist.add_action(
                    Action_Dihedral(), cm, _top,
                    dslist=dslist, *args, **kwd)
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
        for idx, frame in enumerate(_frame_iter_master(traj)):
            arr[idx] = frame.calc_dihedral(int_2darr)

        arr = arr.T
        if dtype == 'ndarray':
            return arr
        else:
            from pytraj.datasetlist import from_dict
            py_dslist = from_dict({'dihedral': arr})
            return _get_data_from_dtype(py_dslist, dtype)


def calc_mindist(traj=None,
                 command="",
                 top=None,
                 dtype='ndarray', *args, **kwd):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.mindist(traj, '@CA @H')
    '''
    from pytraj.actions.CpptrajActions import Action_NativeContacts
    from pytraj.utils.convert import array2d_to_cpptraj_maskgroup
    act = Action_NativeContacts()
    dslist = CpptrajDatasetList()

    if not isinstance(command, string_types):
        command = array2d_to_cpptraj_maskgroup(command)
    _command = "mindist " + command
    _top = _get_top(traj, top)
    act(_command, traj, top=_top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)[-1]


def calc_watershell(traj=None, solute_mask=None,
                    solvent_mask=':WAT',
                    lower=3.4, upper=5.0,
                    image=True,
                    dtype='dataset', top=None):
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
    >>> pt.watershell(traj, solute_mask='!:WAT')
    >>> pt.watershell(traj, solute_mask='!:WAT', lower=5.0, upper=10.)
    """
    from pytraj.actions.CpptrajActions import Action_Watershell
    _top = _get_top(traj, top)
    _solutemask = solute_mask if solute_mask is not None else ''

    if _solutemask in [None, '']:
        raise ValueError('must provide solute mask')

    _solventmask = solvent_mask if solvent_mask is not None else ''
    _noimage= 'noimage' if not image else ''

    _lower = 'lower ' + str(lower)
    _upper = 'upper ' + str(upper)
    command = ' '.join((_solutemask, _lower, _upper, _noimage, _solventmask))

    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    if not 'out' in command:
        # current Watershell action require specifying output
        command += ' out tmp.tmp'
    dslist = CpptrajDatasetList()
    Action_Watershell()(command, traj, _top, dslist=dslist)
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_radial(traj=None, command="", top=Topology()):
    '''Action_Radial require calling Print() to get output. We make change here'''
    act = adict['radial']
    # add `radial` keyword to command (need to check `why`?)
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    command = 'radial ' + command
    dslist = CpptrajDatasetList()
    if not top.is_empty():
        act(command, traj, top, dslist=dslist)
    else:
        act(command, traj, dslist=dslist)

    # dump data to dslist.
    act.print_output()
    return dslist


def calc_matrix(traj=None,
                command="",
                top=None,
                dtype='ndarray', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Matrix
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    act = Action_Matrix()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype)


def calc_radgyr(traj=None,
                mask="",
                top=None,
                nomax=True,
                dtype='ndarray', *args, **kwd):
    '''calc radgyr

    Examples
    --------
    >>> pt.radgyr(traj, '@CA')
    >>> pt.radgyr(traj, '!:WAT', nomax=False)
    '''

    from pytraj.actions.CpptrajActions import Action_Radgyr

    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    _nomax = 'nomax' if nomax else ""
    command = " ".join((mask, _nomax))

    act = Action_Radgyr()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_molsurf(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    '''calc molsurf

    Examples
    --------
    >>> pt.molsurf(traj, '@CA')
    >>> pt.molsurf(traj, '!:WAT')
    '''
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)
    command = mask

    from pytraj.actions.CpptrajActions import Action_Molsurf
    act = Action_Molsurf()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_distrmsd(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = mask

    from pytraj.actions.CpptrajActions import Action_DistRmsd
    act = Action_DistRmsd()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_volume(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = mask

    from pytraj.actions.CpptrajActions import Action_Volume
    act = Action_Volume()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_multivector(traj=None,
                     command="",
                     top=None,
                     dtype='ndarray', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_MultiVector
    act = Action_MultiVector()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_volmap(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = mask

    from pytraj.actions.CpptrajActions import Action_Volmap
    act = Action_Volmap()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def calc_linear_interaction_energy(traj=None,
                                   mask="",
                                   top=None,
                                   dtype='dataset', *args, **kwd):
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = mask
    from pytraj.actions.CpptrajActions import Action_LIE
    act = Action_LIE()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)

# alias
calc_LIE = calc_linear_interaction_energy


def calc_rdf(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    from pytraj.actions.CpptrajActions import Action_Radial
    act = Action_Radial()
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)

    command = "pytraj_tmp_output.agr " + command
    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype)


def calc_jcoupling(traj=None,
                   mask="",
                   top=None,
                   kfile=None,
                   dtype='dataset', *args, **kwd):
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
    """
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)
    command = mask

    from pytraj.actions.CpptrajActions import Action_Jcoupling
    act = Action_Jcoupling()
    # add `radial` keyword to command (need to check `why`?)
    dslist = CpptrajDatasetList()
    _top = _get_top(traj, top)

    if kfile is not None:
        command += " kfile %s" % kfile
    act(command, traj, dslist=dslist, top=_top, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype)


def do_translation(traj=None, command="", top=None):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> # load to mutable trajectory by `load` method
    >>> traj = pt.load('traj.nc', 'myparm.parm7')
    >>> pt.translate(traj, '@CA x 120.')
    '''
    from pytraj.actions.CpptrajActions import Action_Translate

    _noaction_with_TrajectoryIterator(traj)

    _top = _get_top(traj, top)

    if is_array(command):
        x, y, z = command
        _x = "x " + str(x)
        _y = "y " + str(y)
        _z = "z " + str(z)
        _command = " ".join((_x, _y, _z))
    else:
        _command = command
    Action_Translate()(_command, traj, top=_top)


translate = do_translation


def do_scaling(traj=None, command="", top=None):
    from pytraj.actions.CpptrajActions import Action_Scale
    _noaction_with_TrajectoryIterator(traj)
    _top = _get_top(traj, top)
    Action_Scale()(command, traj, top=_top)

scale = do_scaling

def do_rotation(traj=None, command="", top=None):
    from pytraj.actions.CpptrajActions import Action_Rotate
    _top = _get_top(traj, top)
    _noaction_with_TrajectoryIterator(traj)
    Action_Rotate()(command, traj, top=_top)


rotate = do_rotation


def do_autoimage(traj=None, command="", top=None):
    _noaction_with_TrajectoryIterator(traj)
    _top = _get_top(traj, top)
    from pytraj.actions.CpptrajActions import Action_AutoImage
    Action_AutoImage()(command, traj, top=_top)

autoimage = do_autoimage


def get_average_frame(traj=None, command="", top=None):
    from pytraj.actions.CpptrajActions import Action_Average
    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)

    # add "crdset s1" to trick cpptraj dumpt coords to DatSetList
    command += " crdset s1"

    act = Action_Average()
    act(command, traj, _top, dslist=dslist)

    # need to call this method so cpptraj will write
    act.print_output()

    return dslist[0].get_frame()


def randomize_ions(traj=None, command="", top=None):
    """randomize_ions for given Frame with Topology

    Parameters
    ----------
    traj : Trajectory-like or a Frame
    top : Topology, optional (only needed if ``traj`` does not have it)

    Notes
    -----
    ``traj`` must be mutable since this method inplace update coordinate

    """
    from pytraj.actions.CpptrajActions import Action_RandomizeIons
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    _noaction_with_TrajectoryIterator(traj)
    act = Action_RandomizeIons()
    act(command, traj, top)


def clustering_dataset(array_like, command=''):
    '''
    Returns
    -------
    cluster index for each data point

    Examples
    --------
    >>> import pytraj as pt
    >>> import numpy as np 
    >>> array_like = np.random.randint(0, 10, 1000)
    >>> data = pt.clustering_dataset(array_like, 'clusters 10 epsilon 3.0')
    >>> print(data)
    '''
    import numpy as np
    from pytraj.analyses.CpptrajAnalyses import Analysis_Clustering

    dslist = CpptrajDatasetList()
    dslist.add_set('double', '__array_like')
    dslist[0].resize(len(array_like))
    dslist[0].values[:] = array_like
    act = Analysis_Clustering()
    command = 'data __array_like ' + command
    act(command, dslist=dslist)

    return np.array(dslist[-1])


def do_clustering(traj=None,
                  command="",
                  top=None,
                  dtype='dataset',
                  dslist=None,
                  dflist=None):
    """
    Parameters
    ---------
    traj : Trajectory-like | list of Trajectory-like | frame or chunk iterator
    command : cpptraj command
    top : Topology, optional
    dslist : CpptrajDatasetList, optional
    dflist : DataFileList, optional

    Notes:
    Supported algorithms: kmeans, hieragglo, and dbscan.

    Examples
    --------
        do_clustering(traj, "kmeans clusters 50 @CA")

    Returns
    -------
    CpptrajDatasetList object

    """

    _top = _get_top(traj, top)
    ana = analdict['clustering']
    # need to creat `dslist` here so that every time `do_clustering` is called,
    # we will get a fresh one (or will get segfault)
    if dslist is None:
        dslist = CpptrajDatasetList()
    else:
        dslist = dslist

    if traj is not None:
        dslist.add_set("coords", "__pytraj_cluster")
        #dslist[-1].top = _top
        dslist[0].top = _top
        for frame in traj:
            # dslist[-1].add_frame(frame)
            dslist[0].add_frame(frame)
        command += " crdset __pytraj_cluster"
    else:
        pass
    ana(command, _top, dslist, dflist)
    # remove frames in dslist to save memory
    dslist.remove_set(dslist['__pytraj_cluster'])
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_multidihedral(traj=None,
                       dhtypes=None,
                       resrange=None,
                       define_new_type=None,
                       range360=False,
                       dtype='dataset',
                       top=None, *args, **kwd):
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
    >>> pt.multidihedral(traj)
    >>> pt.multidihedral(traj, 'phi psi')
    >>> pt.multidihedral(traj, resrange=range(8))
    >>> pt.multidihedral(traj, range360=True)
    >>> pt.multidihedral(traj, resrange='1,3,5')
    >>> pt.multidihedral(traj, dhtypes='phi psi')
    >>> pt.multidihedral(traj, dhtypes='phi psi', resrange='3-7')
    >>> pt.multidihedral(traj, dhtypes='phi psi', resrange=[3, 4, 8])
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

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act = adict['multidihedral']
    act(_command, traj, _top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_atomicfluct(traj=None,
                     mask="",
                     top=None,
                     dtype='dataset', *args, **kwd):
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = mask

    _top = _get_top(traj, top)

    dslist = CpptrajDatasetList()
    act = adict['atomicfluct']
    act(command, traj, top=_top, dslist=dslist, *args, **kwd)
    # tag: print_output()
    act.print_output()  # need to have this. check cpptraj's code
    return _get_data_from_dtype(dslist, dtype=dtype)


@noparallel
def calc_bfactors(traj=None,
                  mask="",
                  byres=True,
                  top=None,
                  dtype='ndarray', *args, **kwd):
    """
    Returns
    -------
    if dtype is 'ndarray' (default), return a numpy array
    with shape=(n_atoms/n_residues, 2) ([atom_or_residue_idx, value])

    Examples
    --------
    >>> import pytraj as pt
    >>> pt.calc_bfactors(traj, byres=True)
    """
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    byres_text = "byres" if byres else ""

    _command = " ".join((mask, byres_text, "bfactor"))
    return calc_atomicfluct(traj=traj,
                            mask=_command,
                            top=top,
                            dtype=dtype, *args, **kwd)


def calc_vector(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """perform dihedral search
    Parameters
    ----------
    command : str, cpptraj command 
    traj : Trajectory-like object
    *arg and **kwd: additional arguments

    Returns
    -------
    DataSet_Vector object

    Examples
    ------
    >>> import pytraj.common_actions as pyca
    >>> pyca.calc_vector(traj, "@CA @CB").tolist()
    >>> pyca.calc_vector(traj, "", traj).tolist()
    >>> pyca.calc_vector(traj, "principal z").to_ndarray()
    >>> pyca.calc_vector(traj, "principal x").to_ndarray()
    >>> pyca.calc_vector(traj, "ucellx").tolist()
    >>> pyca.calc_vector(traj, "boxcenter").tolist()
    >>> pyca.calc_vector(traj, "box").tolist()
    """
    from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
    from pytraj.actions.CpptrajActions import Action_Vector
    from pytraj.core.ActionList import ActionList

    dslist = CpptrajDatasetList()
    _top = _get_top(traj, top)
    list_of_commands = _get_list_of_commands(mask)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        actlist.add_action(act, command, _top, dslist=dslist, *args, **kwd)
    actlist.do_actions(traj)

    if dtype == 'vector':
        return dslist[-1]
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def _calc_vector_center(traj=None,
                        command="",
                        top=None,
                        mass=False,
                        dtype='ndarray'):
    from pytraj.actions.CpptrajActions import Action_Vector
    _top = _get_top(traj, top)

    dslist = CpptrajDatasetList()
    dslist.set_py_free_mem(False)  # need this to avoid segmentation fault
    act = Action_Vector()
    command = "center " + command

    if mass:
        command += " mass"

    act.read_input(command=command, top=_top, dslist=dslist)
    act.process(_top)

    for frame in _frame_iter_master(traj):
        # set Frame masses
        if mass:
            frame.set_frame_mass(_top)
        act.do_action(frame)
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_center_of_mass(traj=None,
                        mask='',
                        top=None,
                        dtype='ndarray', *args, **kwd):
    return _calc_vector_center(traj=traj,
                               command=mask,
                               top=top,
                               mass=True,
                               dtype=dtype)


calc_COM = calc_center_of_mass


def calc_center_of_geometry(traj=None, command="", top=None, dtype='ndarray'):
    _top = _get_top(traj, top)

    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)

    atom_mask_obj = _top(command)
    dslist = CpptrajDatasetList()
    dslist.add_set("vector")

    for frame in _frame_iter_master(traj):
        dslist[0].append(frame.center_of_geometry(atom_mask_obj))
    return _get_data_from_dtype(dslist, dtype=dtype)


calc_COG = calc_center_of_geometry


def calc_pairwise_rmsd(traj=None,
                       mask="",
                       metric='rms',
                       top=None,
                       dtype='ndarray',
                       mat_type='full', *args, **kwd):
    """return 2D numpy array

    Parameters
    ----------
    traj : Trajectory-like, iterable object
    mask : mask (default=all atom) + extra command
        See `Notes` below for further info
    metric : {'rms', 'dme', 'srmsd', 'nofit'}
    top : Topology, optional, default=None
    dtype: ndarray
    mat_type : return matrix type, default 'full'
    *args, **kwd: optional (for advanced user)

    Examples
    --------
    >>> pt.pairwise_rmsd(traj(0, 1000, mask='@CA'))

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
    This calculation is memory consumming. It's better to use **pytraj.TrajectoryIterator**
    (**pytraj.iterload(...)**)

    It's better to use `**pytraj.pairwise_rmsd(traj(mask='@CA'))**` than
    `**pytraj.pairwise_rmsd(traj, mask='@CA')**`.

    Install **libcpptraj** with openmp to benifit from parallel
    """
    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    command = ' '.join((mask, metric))

    from pytraj.analyses.CpptrajAnalyses import Analysis_Rms2d
    act = Analysis_Rms2d()

    dslist = CpptrajDatasetList()
    dslist.add_set("coords", "mylovelypairwisermsd")
    _top = _get_top(traj, top)
    dslist[0].top = _top
    # need to set "rmsout" to trick cpptraj not giving error
    # need " " (space) before crdset too
    command = command + " crdset mylovelypairwisermsd rmsout mycrazyoutput"

    # upload Frame to crdset
    for frame in _frame_iter_master(traj):
        dslist[0].append(frame)

    act(command, _top, dslist=dslist, *args, **kwd)
    # remove dataset coords to free memory
    dslist.remove_set(dslist[0])

    if dtype == 'ndarray':
        return _get_matrix_from_dataset(dslist[0], mat_type)
    else:
        return _get_data_from_dtype(dslist, dtype)


def calc_density(traj=None,
                 command="",
                 top=None,
                 dtype='ndarray', *args, **kwd):
    # NOTE: trick cpptraj to write to file first and the reload

    with goto_temp_folder():

        def _calc_density(traj, command, *args, **kwd):
            # TODO: update this method if cpptraj save data to
            # CpptrajDatasetList
            from pytraj.actions.CpptrajActions import Action_Density

            _top = _get_top(traj, top)
            dflist = DataFileList()

            tmp_filename = "tmp_pytraj_out.txt"
            command = "out " + tmp_filename + " " + command
            act = Action_Density()
            # with goto_temp_folder():
            act(command, traj, top=_top, dflist=dflist)
            act.print_output()
            dflist.write_all_datafiles()
            absolute_path_tmp = os.path.abspath(tmp_filename)
            return absolute_path_tmp

        dslist = CpptrajDatasetList()
        fname = _calc_density(traj, command, *args, **kwd)
        dslist.read_data(fname)
        return _get_data_from_dtype(dslist, dtype)


def calc_temperatures(traj=None, command="", top=None, dtype='ndarray'):
    """return 1D python array of temperatures (from velocity) in traj
    if `frame` keyword is specified cpptraj/pytraj will take existing T

    Default = array of 0.0
    """
    from array import array as pyarray
    _top = _get_top(traj, top)
    dslist = calculate('temperature', traj, command, _top)
    return _get_data_from_dtype(dslist, dtype)


def rmsd_perres(traj=None,
                ref=0,
                mask="",
                mass=False,
                top=None,
                range=None,
                perresmask=None,
                dtype='dataset', *args, **kwd):
    """
    Perform rmsfit calculation with `mask`, then calculate nofit rms for residues
    in `range` with given `perresmask`
    """
    from pytraj.utils.convert import array_to_cpptraj_range
    if range is not None:
        if isinstance(range, string_types):
            _range = 'range %s ' % range
        else:
            raise ValueError("range must be a string")
    else:
        _range = ''
    _perresmask = perresmask if perresmask is not None else ''
    cm = " ".join((mask, 'perres', _range, _perresmask))
    return calc_rmsd(traj=traj,
                     ref=ref,
                     mask=cm,
                     nofit=False,
                     mass=mass,
                     top=top,
                     dtype=dtype, *args, **kwd)


def calc_rmsd(traj=None,
              ref=0,
              mask="",
              nofit=False,
              mass=False,
              top=None,
              dtype='ndarray', *args, **kwd):
    """calculate rmsd

    Parameters
    ----------
    traj : Trajectory-like | List of trajectories | Trajectory or frame_iter
    ref : {Frame, int}, default=0
        Reference frame or index.
    mask : str or 1D array-like of string or 1D or 2D array-like
        Atom mask/indices
    top : {Topology, str}, default None, optional
    dtype : return data type, default='ndarray'

    Examples
    --------
    >>> import pytraj as pt
    >>> pt.rmsd(traj, ref=-3) # ref=traj[-3]
    >>> pt.rmsd(traj, mask=['@CA', '@C', ':3-18@CA'], dtype='dataset')
    >>> pt.rmsd(traj, ref=traj[0], mask=':3-13')

    """
    from pytraj.utils import is_int
    from array import array as pyarray
    from pytraj.datasets import DatasetDouble
    from pytraj.actions.CpptrajActions import Action_Rmsd
    from pytraj.core.ActionList import ActionList
    import numpy as np

    _nofit = ' nofit ' if nofit else ''
    _mass = ' mass ' if mass else ''
    opt = _nofit + _mass

    if isinstance(mask, string_types):
        command = [mask, ]
    else:
        try:
            cmd = np.asarray(mask)
        except ValueError as e:
            raise ValueError("don't mix different types")
        dname = cmd.dtype.name
        if 'str' in dname:
            command = cmd
        elif 'int' in dname or 'object' in dname:
            if cmd.ndim == 1 and 'object' not in dname:
                command = [to_cpptraj_atommask(mask), ]
            elif cmd.ndim == 2 or 'object' in dname:
                command = [to_cpptraj_atommask(x) for x in mask]
            else:
                raise ValueError("only support array with ndim=1,2")
        else:
            raise ValueError("not supported")

    _top = _get_top(traj, top)

    ref = _get_reference_from_traj(traj, ref)

    alist = ActionList()
    dslist = CpptrajDatasetList()

    for cm in command:
        _cm = cm + opt
        alist.add_action(Action_Rmsd(), _cm, top=_top, dslist=dslist)

    alist.do_actions(ref)
    alist.do_actions(traj)

    if dtype == 'pyarray':
        return pyarray('d', dslist[0].data)[1:]
    else:
        from pytraj.datasetlist import DatasetList
        dnew = DatasetList(dslist)
        for d in dnew:
            d.values = d.values[1:]
        return _get_data_from_dtype(dnew, dtype=dtype)

# alias for `calc_rmsd`
rmsd = calc_rmsd


def calc_rmsd_with_rotation_matrices(
    traj=None,
    mask="",
    ref=None,
    top=None,
    dtype='dataset', *args, **kwd):

    ref = _get_reference_from_traj(traj, ref)

    if not isinstance(mask, string_types):
        # [1, 3, 5] to "@1,3,5
        mask = to_cpptraj_atommask(mask)

    command = mask

    if not isinstance(command, string_types):
        raise ValueError("only support string mask/command in mode=cpptraj")
    if dtype in ['ndarray', 'pyarray']:
        raise ValueError("does not support ndarray/pyarray here")

    ref = _get_reference_from_traj(traj, ref)

    _top = _get_top(traj, top)
    from pytraj.actions.CpptrajActions import Action_Rmsd
    act = Action_Rmsd()
    dslist = CpptrajDatasetList()
    act(command + " savematrices", [ref, traj], top=_top, dslist=dslist)

    # skip reference frame
    from pytraj.datasetlist import DatasetList
    dslist = DatasetList(dslist)
    dslist[0].values = dslist[0].values[1:]
    dslist[1].values = dslist[1].values[1:]
    return _get_data_from_dtype(dslist, dtype=dtype)


def align_principal_axis(traj=None, mask="*", top=None):
    # TODO : does not match with cpptraj output
    # rmsd_nofit ~ 0.5 for md1_prod.Tc5b.x, 1st frame
    """
    Notes
    -----
    apply for mutatble traj (Trajectory, Frame)
    """
    from pytraj.actions.CpptrajActions import Action_Principal
    _noaction_with_TrajectoryIterator(traj)

    command = mask

    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    _top = _get_top(traj, top)
    act = Action_Principal()
    command += " dorotation"
    act(command, traj, top=_top)


def principal_axes(traj=None,
                   mask='*',
                   top=None,
                   dorotation=False,
                   mass=True,
                   dtype='dataset', *args, **kwd):
    """
    Returns
    -------
    pytraj.DatasetList of matrix with shape=(n_frames, 3, 3) and vector with shape=(n_frames, 3)
    if `dorotation`, the system will be aligned along principal axes
    (apply for mutable system)
    """
    from pytraj.actions.CpptrajActions import Action_Principal
    act = Action_Principal()
    command = mask

    _dorotation = 'dorotation' if dorotation else ''
    _mass = 'mass' if mass else ''

    if 'name' not in command:
        command += ' name pout'

    command = ' '.join((command, _dorotation, _mass))
    print(command)

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()
    act(command, traj, _top, dslist=dslist, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)


def atomiccorr(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """
    """
    from pytraj.actions.CpptrajActions import Action_AtomicCorr
    _top = _get_top(traj, top)

    if not isinstance(mask, string_types):
        # [1, 3, 5] to "@1,3,5
        mask = to_cpptraj_atommask(mask)

    command = mask

    dslist = CpptrajDatasetList()
    act = adict['atomiccorr']
    act("out mytempfile.out " + command, traj,
        top=_top,
        dslist=dslist, *args, **kwd)
    act.print_output()
    return _get_data_from_dtype(dslist, dtype=dtype)


def _closest_iter(act, traj):
    '''
    Parameters
    ----------
    act : Action object
    traj : Trajectory-like
    '''

    for frame in _frame_iter_master(traj):
        new_frame = Frame()
        new_frame.py_free_mem = False  # cpptraj will do
        act.do_action(frame, new_frame)
        yield new_frame


def closest(traj=None, mask='*', n_solvents=0, restype='trajectory', top=None):
    """return either a new Trajectory or a frame iterator. Keep only ``n_solvents`` closest to mask

    Parameters
    ----------
    traj : Trajectory-like | list of Trajectory-like/frames | frame iterator | chunk iterator
    mask: str, default '*' (all solute atoms)
    restype : str, {'trajectory', 'dataset', 'iterator'}, default 'trajectory'
        if restype == 'trajectory', return a new ``pytraj.Trajectory``
        if restype == 'dataset': return a tuple (new_traj, datasetlist)
        if restype == 'iterator': return a tuple of (Frame iterator, new Topology),  good for memory saving
        if restype == 'all': return (Trajectory, DatasetList)
    top : Topology-like object, default=None, optional

    Returns
    -------
    out : (check above)

    Examples
    --------
    >>> # obtain new traj, keeping only closest 100 waters 
    >>> # to residues 1 to 13 (index starts from 1) by distance to the first atom of water
    >>> t = pt.closest(traj, mask='@CA', n_solvents=10)

    >>> # only get meta data for frames, solvent without getting new Trajectory
    >>> # (such as Frame number, original solvent molecule number, ...) (from cpptraj manual)
    >>> dslist = pt.closest(traj, n_solvents=100, mask=':1-13', restype='dataset')

    >>> # getting a frame iterator for lazy evaluation
    >>> fiter = pt.closest(traj, n_solvents=20, restype='iterator')
    >>> for frame in fiter: print(frame) 
    """

    from .actions.CpptrajActions import Action_Closest
    from pytraj.Trajectory import Trajectory
    dslist = CpptrajDatasetList()

    if n_solvents == 0:
        raise ValueError('must specify the number of solvents')

    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(command)

    command = str(n_solvents) + ' ' + mask

    dtype = restype

    act = Action_Closest()

    _top = _get_top(traj, top)
    new_top = Topology()
    new_top.py_free_mem = False  # cpptraj will do

    if dtype not in ['trajectory', 'iterator']:
        # trick cpptraj to dump data to CpptrajDatasetList too
        command = command + " closestout tmp_pytraj_closestout.out"

    act.read_input(command, _top, dslist=dslist)
    act.process(_top, new_top)

    fiter = _closest_iter(act, traj)

    if dtype == 'iterator':
        return (fiter, new_top)
    else:
        if dtype in ['trajectory', 'all']:
            fa = Trajectory()
            fa.top = new_top.copy()
            for new_frame in fiter:
                fa.append(new_frame.copy())
            if dtype == 'trajectory':
                return fa
            elif dtype == 'all':
                return fa, dslist
        else:
            for new_frame in fiter:
                # just let cpptraj dump data to DatasetList
                pass
            new_dslist = _get_data_from_dtype(dslist, dtype=dtype)
            return new_dslist


def native_contacts(traj=None,
                    mask="",
                    dtype='dataset',
                    ref=0,
                    distance=7.0,
                    noimage=False,
                    include_solvent=False,
                    byres=False,
                    top=None,
                    *args, **kwd):
    """
    Examples
    --------
    >>> # use 1st frame as reference, don't need specify ref Frame
    >>> pt.native_contacts(traj)

    >>> # explicitly specify reference, specify distance cutoff
    >>> pt.native_contacts(traj, ref=ref, distance=8.0)

    Notes
    -----
    if `ref` is None: first number in result corresponds to reference
    Not assert to cpptraj's output yet
    """
    from .actions.CpptrajActions import Action_NativeContacts
    act = Action_NativeContacts()
    dslist = CpptrajDatasetList()

    if not isinstance(mask, string_types):
        # [1, 3, 5] to "@1,3,5
        mask = to_cpptraj_atommask(mask)

    command = mask

    ref = _get_reference_from_traj(traj, ref)

    _distance = str(distance)
    _noimage = "noimage" if noimage else ""
    _includesolvent = "includesolvent" if include_solvent else ""
    _byres = "byresidue" if byres else ""

    _command = " ".join((command, _distance, _noimage, _includesolvent, _byres
                         ))

    _top = _get_top(traj, top)
    act(_command, [ref, traj], top=_top, dslist=dslist, *args, **kwd)

    from pytraj.datasetlist import DatasetList
    dslist = DatasetList(dslist)
    for d in dslist:
        # exclude ref frame
        d.values = d.values[1:]
    return _get_data_from_dtype(dslist, dtype=dtype)


def calc_grid(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """
    """
    # TODO: doc, rename method, move to seperate module?
    from .actions.CpptrajActions import Action_Grid
    if not isinstance(command, string_types):
        command = to_cpptraj_atommask(command)
    act = Action_Grid()
    dslist = CpptrajDatasetList()

    # cpptraj require output
    command = "tmp_pytraj_grid_output.txt " + command
    _top = _get_top(traj, top)
    act(command, traj, dslist=dslist, top=_top, *args, **kwd)
    return _get_data_from_dtype(dslist, dtype=dtype)


calc_grid.__doc__ = _long_manual.__grid__


def check_structure(traj=None, command="", top=None, *args, **kwd):
    """
    Examples
    --------
    >>> check_structure(traj[0], top=traj.top)
    """
    from .actions.CpptrajActions import Action_CheckStructure
    act = Action_CheckStructure()

    # cpptraj require output
    _top = _get_top(traj, top)
    act(command, traj, top=_top, *args, **kwd)


def timecorr(vec0, vec1,
             order=2,
             timestep=1.,
             tcorr=10000.,
             norm=False,
             dtype='ndarray'):
    """TODO: doc. not yet assert to cpptraj's output

    Parameters
    ----------
    vec0 : 2D array-like, shape=(n_frames, 3)
    vec1 : 2D array-like, shape=(n_frames, 3)
    order : int, default 2
    timestep : float, default 1.
    tcorr : float, default 10000.
    norm : bool, default False
    """
    from pytraj.datasets import DataSetList as CDSL, DatasetVector
    from pytraj.math import Vec3
    import numpy as np
    act = analdict['timecorr']

    cdslist = CDSL()

    cdslist.add_set("vector", "_vec0")
    cdslist.add_set("vector", "_vec1")
    cdslist[0].from_array_like(np.asarray(vec0).astype('f8'))
    cdslist[1].from_array_like(np.asarray(vec1).astype('f8'))

    _order = "order " + str(order)
    _tstep = "tstep " + str(timestep)
    _tcorr = "tcorr " + str(tcorr)
    _norm = "norm" if norm else ""
    command = " ".join(
        ('vec1 _vec0 vec2 _vec1', _order, _tstep, _tcorr, _norm))
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
    from pytraj.datasets import DataSetList as CDSL
    from pytraj.analyses.CpptrajAnalyses import Analysis_CrankShaft
    import numpy as np

    cdslist = CDSL()
    cdslist.add_set("double", "d0")
    cdslist.add_set("double", "d1")

    cdslist[0].from_array_like(np.asarray(data0))
    cdslist[1].from_array_like(np.asarray(data1))

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
    from pytraj.datasets import DataSetList as CDSL
    import numpy as np

    cdslist = CDSL()
    cdslist.add_set("double", "d0")
    cdslist.add_set("double", "d1")

    cdslist[0].from_array_like(np.asarray(data0))
    cdslist[1].from_array_like(np.asarray(data1))

    act = analdict['corr']
    act("d0 d1 out _tmp.out", dslist=cdslist)
    return _get_data_from_dtype(cdslist[2:], dtype=dtype)


def auto_correlation_function(data, dtype='ndarray', covar=True):
    """
    Notes
    -----
    Same as `autocorr` in cpptraj
    """
    from pytraj.datasets import DataSetList as CDSL
    import numpy as np

    _nocovar = " " if covar else " nocovar"

    cdslist = CDSL()
    cdslist.add_set("double", "d0")

    cdslist[0].from_array_like(np.asarray(data))

    act = analdict['autocorr']
    command = "d0 out _tmp.out" + _nocovar
    act(command, dslist=cdslist)
    return _get_data_from_dtype(cdslist[1:], dtype=dtype)


def lifetime(data, command="", dtype='ndarray', *args, **kwd):
    """
    Notes
    -----
    Same as `autocorr` in cpptraj
    """
    from pytraj.datasets import DataSetList as CDSL
    from pytraj.analyses.CpptrajAnalyses import Analysis_Lifetime
    import numpy as np

    cdslist = CDSL()
    if 'int' in data.dtype.name:
        cdslist.add_set("integer", "d0")
    else:
        cdslist.add_set("double", "d0")

    cdslist[0].from_array_like(np.asarray(data))

    act = Analysis_Lifetime()
    command = " ".join((command, "d0"))
    act(command, dslist=cdslist)
    return _get_data_from_dtype(cdslist[1:], dtype=dtype)


def search_neighbors(traj=None,
                      mask='',
                      cutoff='',
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
    >>> pt.search_neighbors(traj, mask=':5', cutoff="<5.0") # around residue 5 with 5.0 cutoff
    >>> pt.search_neighbors(traj, [3, 7, 8], cutoff=">10.0") # around atom 3, 7, 8, larger than 10.0 Angstrom
    """
    from pytraj.datasetlist import DatasetList
    import numpy as np

    if not isinstance(mask, string_types):
        mask = to_cpptraj_atommask(mask)

    if '>' in cutoff:
        cutoff = '>:' + cutoff.split('>')[-1]
    elif '<' in cutoff:
        cutoff = '<:' + cutoff.split('<')[-1]
    else:
        raise ValueError('must correctly specify cutoff: using > or <')

    mask = " ".join((mask, cutoff))

    dslist = DatasetList()

    _top = _get_top(traj, top)
    for idx, frame in enumerate(_frame_iter_master(traj)):
        _top.set_reference_frame(frame)
        dslist.append({str(idx): np.asarray(_top.select(mask))})
    return _get_data_from_dtype(dslist, dtype)


def pucker(traj=None,
           pucker_mask=("C1'", "C2'", "C3'", "C4'", "O4'"),
           resrange=None,
           top=None,
           dtype='dataset',
           range360=False,
           method='altona',
           use_com=True,
           amplitude=True,
           offset=None, *args, **kwd):
    """Note: not validate yet

    """
    from pytraj.datasetlist import DatasetList
    from pytraj.datasets import DataSetList as CDL
    from pytraj.actions.CpptrajActions import Action_Pucker
    from pytraj.compat import range

    _top = _get_top(traj, top)
    if not resrange:
        resrange = range(_top.n_residues)

    _range360 = "range360" if range360 else ""
    geom = "geom" if not use_com else ""
    amp = "amplitude" if amplitude else ""
    _offset = "offset " + str(offset) if offset else ""

    cdslist = CDL()
    for res in resrange:
        act = Action_Pucker()
        command = " ".join((":" + str(res + 1) + '@' + x for x in pucker_mask))
        name = "pucker_res" + str(res + 1)
        command = " ".join((name, command, _range360, method, geom, _offset))
        act(command, traj, top=_top, dslist=cdslist, *args, **kwd)
    return _get_data_from_dtype(cdslist, dtype)


def center(traj=None, mask="", top=None):
    """
    Examples
    --------
    >>> pt.center(traj) # all atoms, center to box center (x/2, y/2, z/2)
    >>> pt.center(traj, '@CA origin') # center at origin, use @CA
    >>> pt.center(traj, 'mass') # center to box center, use mass weighted.
    >>> pt.center(traj, ':1 mass') # residue 1, use mass weighted.
    """
    _noaction_with_TrajectoryIterator(traj)
    _top = _get_top(traj, top)
    from pytraj.actions.CpptrajActions import Action_Center
    act = Action_Center()
    act(mask, traj, top=_top)


def rotate_dihedral(traj=None, mask="", top=None):
    # change to pt.rotate_dihedral(traj, res=0, 
    #              mask=("O4'", "C1'", "N9", "C4"), deg=120)?
    """
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.rotate_dihedral(traj, "3:phi:120") # rotate phi of res 3 to 120 deg
    >>> pt.rotate_dihedral(traj, "1:O4':C1':N9:C4:120") # rotate dihedral with given mask

    Notes
    -----
    Syntax and method's name might be changed
    """
    _noaction_with_TrajectoryIterator(traj)
    _top = _get_top(traj, top)

    if "custom:" in mask:
        command = mask
    else:
        command = "custom:" + mask

    from pytraj.actions.CpptrajActions import Action_MakeStructure
    act = Action_MakeStructure()

    act(command, traj, top=_top)


def _rotate_dih(traj, resid='1', dihtype=None, deg=0, top=None):
    '''
    Examples
    >>> pt._rotate_dih(traj, resid='1', dihtype='delta')
    '''
    _top = _get_top(traj, top)

    if not isinstance(resid, string_types):
        resid = str(resid + 1)
    deg = str(deg)

    command = ':'.join((dihtype, resid, dihtype, deg))
    make_structure(traj, command, top=_top)


set_dihedral = _rotate_dih


def make_structure(traj=None, mask="", top=None):
    from pytraj.actions.CpptrajActions import Action_MakeStructure
    _noaction_with_TrajectoryIterator(traj)
    _top = _get_top(traj, top)

    command = mask
    act = Action_MakeStructure()
    act(command, traj, top=_top)
