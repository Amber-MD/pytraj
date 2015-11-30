from __future__ import absolute_import
# TODO: rename this file
from functools import wraps

# do not import anything else here.
from pytraj.externals.six import string_types, integer_types
from pytraj.utils.convert import array_to_cpptraj_atommask


def _load_Topology(filename):
    from pytraj import Topology, ParmFile
    top = Topology()
    parm = ParmFile()
    parm.readparm(filename, top)
    return top


def _get_topology(traj, top):
    '''

    >>> import pytraj as pt
    >>> from pytraj import Topology
    >>> traj = pt.datafiles.load_rna()
    >>> top_filename = traj.top.filename
    >>> isinstance(_get_topology(traj, None), Topology)
    True
    >>> isinstance(_get_topology(traj, top_filename), Topology)
    True
    >>> top = traj.top
    >>> isinstance(_get_topology(traj, top), Topology)
    True
    >>> # find Topology in a list
    >>> isinstance(_get_topology([traj, traj], None), Topology)
    True
    '''
    if isinstance(top, string_types):
        # if provide a filename, load to Topology
        _top = _load_Topology(top)
    elif top is None:
        # if user does not provide Topology, try to find it in traj
        if hasattr(traj, 'top'):
            _top = traj.top
        else:
            # list, tuple of traj objects
            try:
                for tmp in traj:
                    if hasattr(tmp, 'top'):
                        _top = tmp.top
                        break
            except TypeError:
                _top = None
    else:
        _top = top
    return _top


def _get_data_from_dtype(d0, dtype='dataset'):
    from pytraj.datasetlist import DatasetList as DSL

    if dtype is None or dtype == 'dataset':
        if hasattr(d0, 'set_own_memory'):
            d0.set_own_memory(False)

    dtype = dtype.lower()
    if dtype == 'dataset':
        return DSL(d0)
    elif dtype == 'ndarray':
        return d0.to_ndarray()
    elif dtype == 'dict':
        return d0.to_dict()
    elif dtype == 'dataframe':
        if hasattr(d0, 'key'):
            d0.key = d0.key.replace(':', '_')
            d0.key = d0.key.replace('-', '_')
        else:
            for _d in d0:
                _d.key = _d.key.replace(':', '_')
                _d.key = _d.key.replace('-', '_')
        return d0.to_dataframe()
    elif dtype == 'cpptraj_dataset':
        return d0
    else:
        raise NotImplementedError()


def _get_list_of_commands(mask_or_commands):
    '''

    Examples
    --------
    >>> _get_list_of_commands('@CA')
    ['@CA']
    >>> _get_list_of_commands(('@CA'))
    ['@CA']
    >>> _get_list_of_commands(100)
    Traceback (most recent call last):
        ...
    ValueError: must be string or list/tuple of strings
    '''
    if isinstance(mask_or_commands, string_types):
        return [mask_or_commands, ]
    elif isinstance(mask_or_commands, (list, tuple)):
        return list(mask_or_commands)
    else:
        raise ValueError("must be string or list/tuple of strings")


def _get_matrix_from_dataset(dset, mat_type='full'):
    '''return full or half matrix

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> dset = pt.matrix.dist(traj, '@CA', dtype='cpptraj_dataset')[0]
    >>> _get_matrix_from_dataset(dset, mat_type='full').shape
    (12, 12)
    >>> _get_matrix_from_dataset(dset, mat_type='half').shape
    (78,)
    '''
    # dset in DatasetMatrixDouble object
    if mat_type == 'full':
        return dset.values
    elif mat_type == 'half':
        return dset._to_cpptraj_sparse_matrix()
    else:
        raise ValueError()


def _get_reference_from_traj(traj, ref):
    '''try best to get reference

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> frame = _get_reference_from_traj(traj, 3)
    >>> isinstance(frame, pt.Frame)
    True
    >>> frame = _get_reference_from_traj(traj, None)
    >>> isinstance(frame, pt.Frame)
    True
    >>> ref = traj[5]
    >>> frame = _get_reference_from_traj(traj, ref)
    '''
    if isinstance(ref, integer_types):
        try:
            return traj[ref]
        except IndexError:
            raise IndexError("{} does not support indexing".format(traj))
    elif ref is None:
        try:
            return traj[0]
        except IndexError:
            raise IndexError(
                "If reference is an integer, {} must support indexing".format(traj))
    elif 'Trajectory' in ref.__class__.__name__:
        assert ref.n_frames == 1, "only support 1-frame Trajectory as reference"
        return ref[0]
    else:
        return ref


def _get_fiterator(traj, frame_indices=None):
    if frame_indices is None:
        return traj
    else:
        return traj.iterframe(frame_indices=frame_indices)


def _get_resrange(resrange):
    '''return resrange as a string

    Examples
    --------
    >>> _get_resrange('1-3')
    'resrange 1-3'
    >>> _get_resrange(0)
    'resrange 1'
    >>> _get_resrange(range(3))
    'resrange 1,2,3'
    >>> _get_resrange([2, 5, 7])
    'resrange 3,6,8'
    >>> _get_resrange(None)
    ''
    '''
    from pytraj.utils import convert, is_int

    if resrange is not None:
        if is_int(resrange):
            resrange = [resrange, ]
        if isinstance(resrange, string_types):
            _resrange = "resrange " + resrange
        else:
            _resrange = convert.array_to_cpptraj_range(resrange)
            _resrange = "resrange " + str(_resrange)
    else:
        _resrange = ""
    return _resrange


class _super_dispatch(object):
    # TODO: more descriptive method name?
    '''apply a series of functions to ``f``'s args and kwd

    - get Topology from a given traj (Trajectory, frame iterator, ...) and top
        _get_topology(traj, top)
    - create frame iterator from traj and frame_indices
        _get_fiterator(traj, frame_indices)
    - create Amber mask from atom index array
        array_to_cpptraj_atommask(mask)
    - convert int ref to Frame ref
    '''

    def __init__(self, has_ref=False):
        self.has_ref = has_ref

    def __call__(self, f):
        @wraps(f)
        def inner(*args, **kwd):
            args = list(args)
            # traj is always 1st argument
            try:
                traj = kwd.get('traj', args[0])
            except IndexError:
                traj = kwd.get('traj')
            frame_indices = kwd.get('frame_indices')
            ref = kwd.get('ref', None)
            if self.has_ref and ref is None:
                try:
                    ref = args[1]
                except IndexError:
                    ref = 0
            if 'ref' in kwd.keys() or self.has_ref:
                # convert to Frame
                ref = _get_reference_from_traj(traj, ref)

            top = kwd.get('top', None)

            if 'mask' in kwd.keys():
                mask = kwd.get('mask')
                has_mask = True
            else:
                # mask is always 2nd argument if there is no ref
                if not self.has_ref:
                    try:
                        mask = args[1]
                        has_mask = True
                    except IndexError:
                        mask = '*'
                        has_mask = False
                else:
                    try:
                        mask = args[2]
                        has_mask = True
                    except IndexError:
                        mask = '*'
                        has_mask = False

            # overwrite
            kwd['top'] = _get_topology(traj, top)
            if ref is not None:
                if 'ref' in kwd.keys():
                    kwd['ref'] = _get_reference_from_traj(traj, ref)
                else:
                    try:
                        args[1] = ref
                    except IndexError:
                        args.append(ref)
            if 'traj' in kwd.keys():
                kwd['traj'] = _get_fiterator(traj, frame_indices)
            else:
                args[0] = _get_fiterator(traj, frame_indices)
            if not isinstance(mask, string_types):
                mask = array_to_cpptraj_atommask(mask)
            if 'mask' in kwd.keys():
                kwd['mask'] = mask
            else:
                if has_mask:
                    args[1] = mask
            return f(*args, **kwd)

        inner._is_super_dispatched = True
        return inner


def _get_fi_with_dslist(traj, mask, frame_indices, top, crdname='dataset_coords'):
    from pytraj import Trajectory, TrajectoryIterator
    from pytraj.datasets import CpptrajDatasetList
    from pytraj._shared_methods import iterframe_master

    dslist = CpptrajDatasetList()
    dslist.add_set("coords", crdname)
    # need to set "rmsout" to trick cpptraj not giving error
    # need " " (space) before crdset too

    if isinstance(traj, (Trajectory, TrajectoryIterator)):
        # we do atom stripping here before copying to DatasetCoordsCRD to save memory if
        # loading from TrajectoryIterator
        fi = traj.iterframe(mask=mask, frame_indices=frame_indices)
        command = ''
        # use Topology from fi (could be stripped to save memory)
        dslist[0].top = fi.top
        _top = fi.top
    else:
        # ignore frame_indices
        fi = iterframe_master(traj)
        command = mask
        _top = _get_topology(traj, top)
        dslist[0].top = _top
    for frame in fi:
        dslist[0].append(frame)
    return dslist, _top, command
