"""seperate module, only use stdlib
If want to use external package, import it inside the function

This module stores all useful functions that does not fit to anywhere else.
"""
from __future__ import absolute_import
import sys as _sys
import os
from glob import glob
from itertools import islice, groupby
import functools
from collections import OrderedDict


def _array_to_cpptraj_range(seq):
    # use "i+1" since cpptraj use 1-based index for mask
    return ",".join((str(i + 1) for i in seq))


def array_to_atommask(seq):
    '''
    [1, 3, 4] --> @2,4,5
    '''
    return '@' + _array_to_cpptraj_range(seq)


def array_to_atommask_2_groups(seq):
    '''
    [1, 3] --> @1 @3
    [1, 3, 4] --> @1 @3 @4
    '''
    return ' '.join('@' + str(i+1) for i in seq)


def array_to_residuemask(seq):
    '''[1, 3, 4] --> :2,4,5'''
    return ':' + _array_to_cpptraj_range(seq)

# string_types, PY2, PY3, iteritems were copied from six.py
# see license in $PYTRAJHOME/license/externals/
PY2 = _sys.version_info[0] == 2
PY3 = _sys.version_info[0] == 3

if PY3:
    _iteritems = "items"
    string_types = str
else:
    _iteritems = "iteritems"
    string_types = basestring


def iteritems(d, **kw):
    """Return an iterator over the (key, value) pairs of a dictionary."""
    return iter(getattr(d, _iteritems)(**kw))


try:
    # PY3
    from functools import reduce
except ImportError:
    #
    pass

# this module gathers commonly used functions
# from toolz, stackoverflow, ... and from myself
# should make this independent from pytraj

try:
    import numpy as np
except ImportError:
    np = None


def _dispatch_value(func):
    def inner(data, *args, **kwd):
        if hasattr(data, 'values'):
            _data = data.values
        else:
            _data = data
        return func(_data, *args, **kwd)

    inner.__doc__ = func.__doc__
    return inner


def _not_yet_tested(func):
    @functools.wraps(func)
    def inner(*args, **kwd):
        return func(*args, **kwd)

    msg = "This method is not tested. Use it with your own risk"
    inner.__doc__ = "\n".join((func.__doc__, "\n", msg))
    return inner


@_dispatch_value
def split(data, n_chunks_or_array):
    """split `self.data` to n_chunks

    Notes : require numpy (same as `array_split`)
    """
    return np.array_split(data, n_chunks_or_array)


def chunk_average(self, n_chunk, restype='same'):
    '''average by chunk'''
    import numpy as np
    from pytraj.array import DataArray

    data = np.array(list(map(np.mean, split(self, n_chunk))))
    if restype == 'same' and isinstance(self, DataArray):
        new_array = self.shallow_copy()
        new_array.values = data
        return new_array
    else:
        return data


def moving_average(data, n):
    """moving average

    Notes
    -----
    from `stackoverflow <http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python>`_
    """
    window = np.ones(int(n)) / float(n)
    new_data = np.convolve(data, window, 'same')
    if hasattr(data, 'values'):
        new_array = data.shallow_copy()
        new_array.values = new_data
        return new_array
    else:
        return new_data


def pipe(obj, func, *args, **kwargs):
    """Notes: copied from pandas PR
    https://github.com/ghl3/pandas/blob/groupby-pipe/pandas/tools/util.py
    see license in pytraj/license/

    Apply a function to a obj either by
    passing the obj as the first argument
    to the function or, in the case that
    the func is a tuple, interpret the first
    element of the tuple as a function and
    pass the obj to that function as a keyword
    arguemnt whose key is the value of the
    second element of the tuple
    """
    if isinstance(func, tuple):
        func, target = func
        if target in kwargs:
            msg = '%s is both the pipe target and a keyword argument' % target
            raise ValueError(msg)
        kwargs[target] = obj
        return func(*args, **kwargs)
    else:
        return func(obj, *args, **kwargs)


def _compose2(f, g):
    # copied from pandas
    # see license in pytraj/license/
    """Compose 2 callables"""
    return lambda *args, **kwargs: f(g(*args, **kwargs))


def compose(*funcs):
    """
    Notes: copied from pandas (added pytraj's example)
    see license in pytraj/license/

    Compose 2 or more callables

    Examples
    --------
    >>> import pytraj as pt
    >>> func = pt.tools.compose(pt.calc_radgyr, pt.iterload)
    >>> func("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    """
    assert len(funcs) > 1, 'At least 2 callables must be passed to compose'
    return reduce(_compose2, funcs)


def grep(self, key):
    """grep key

    Examples
    --------
    >>> import pytraj as pt
    >>> dslist = pt.calc_multidihedral(traj) 
    >>> pt.tools.grep(dslist, 'psi') 
    """
    new_self = self.__class__()
    for d in self:
        if key in d.key:
            new_self.append(d)
    return new_self


def flatten(x):
    """Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Notes
    -----
    from: http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    Examples
    --------
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        # if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, string_types):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def n_grams(a, n, asarray=False):
    """n_grams

    Parameters
    ----------
    a : sequence
    n : number of elements
    asarray : bool, default False
        if False: return an iterator
        if True: return a numpy array

    Notes
    -----
    adapted from: http://sahandsaba.com/thirty-python-language-features-and-tricks-you-may-not-know.html
    """

    z = (islice(a, i, None) for i in range(n))
    it = zip(*z)

    if not asarray:
        return it
    else:
        import numpy as np
        return np.array([x for x in it])


def dict_to_ndarray(dict_of_array):
    """convert OrderedDict to numpy array

    Examples
    --------
    >>> import pytraj as pt
    >>> dslist = traj.search_hbonds()
    >>> dict_of_array = dslist.to_dict(use_numpy=True)
    >>> np.all(pt.tools.dict_to_ndarray(dict_of_array) == dslist.values)
    True
    """
    if not isinstance(dict_of_array, OrderedDict):
        raise NotImplementedError("support only OrderedDict")
    from pytraj.externals.six import iteritems

    return np.array([v for _, v in iteritems(dict_of_array)])


def concat_dict(iterables):
    """concat dict
    """
    new_dict = {}
    for i, d in enumerate(iterables):
        if i == 0:
            # make a copy of first dict
            new_dict.update(d)
        else:
            for k, v in iteritems(new_dict):
                new_dict[k] = np.concatenate((new_dict[k], d[k]))
    return new_dict


def merge_coordinates(iterables):
    """merge_coordinates from frames
    """
    return np.vstack((f.xyz.copy() for f in iterables))


def merge_frames(iterables):
    """merge from frames to a single Frame. Order matters.
    """
    from pytraj import Frame
    xyz = np.vstack((f.xyz.copy() for f in iterables))
    frame = Frame()
    frame.append_xyz(xyz)
    return frame


def merge_frame_from_trajs(trajlist):
    """
    Examples
    --------
    >>> from frame in pt.tools.merge_frame_from_trajs((traj0, traj1, traj2)):
    >>>     print(frame)
    """
    from pytraj import Frame

    if not isinstance(trajlist, (list, tuple)):
        raise ValueError('input must be a list or tuple of trajectories')
    for iterables in zip(*trajlist):
        yield merge_frames(iterables)


def rmsd_1darray(a1, a2):
    '''rmsd of a1 and a2
    '''
    import numpy as np
    from math import sqrt
    arr1 = np.asarray(a1)
    arr2 = np.asarray(a2)

    if len(arr1.shape) > 1 or len(arr2.shape) > 1:
        raise ValueError("1D array only")

    if arr1.shape != arr2.shape:
        raise ValueError("must have the same shape")

    tmp = sum((arr1 - arr2) ** 2)
    return sqrt(tmp / arr1.shape[0])


def rmsd(a1, a2, flatten=True):
    """rmsd for two array with the same shape

    Parameters
    ----------
    a1, a2: np.ndarray
    flatten : bool, default True
        if True: always flatten two input arrays

    Notes
    -----
    This method is different from ``pytraj.rmsd``
    """
    import numpy as np
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    if a1.shape != a2.shape and not flatten:
        raise ValueError("must have the same shape")
    return rmsd_1darray(a1.flatten(), a2.flatten())


def mean_and_error(a1, a2):
    """calculate mean and error from two 1D array-like
    """
    import numpy as np
    mean = np.mean

    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    assert len(a1.shape) == len(a2.shape) == 1, "1D array"
    return (mean(a1 + a2) / 2, mean(np.abs(a1 - a2)) / 2)


def get_parmed_info(its_obj, att):
    '''for getting info from parmed.Struture'''
    import numpy as np
    return np.asarray([getattr(atom, att) for atom in its_obj.atoms])


def split_parmed_by_residues(struct, start=0, stop=-1, stride=1):
    '''split `ParmEd`'s Structure into different residue
    '''
    from pytraj.compat import range
    from pytraj._cyutils import get_positive_idx

    _stop = get_positive_idx(stop, len(struct.residues))

    for i in range(start, _stop, stride):
        j = ':' + str(i + 1)
        # example: traj[':3']
        yield struct[j]


def split_traj_by_residues(traj, start=0, stop=-1, stride=1):
    '''return a generator

    Examples
    --------
    >>> g = pt.tools.split_traj_by_residues(traj)
    >>> next(g)
    >>> next(g)
    '''
    from pytraj.compat import range
    from pytraj._cyutils import get_positive_idx

    _stop = get_positive_idx(stop, traj.top.n_residues)

    for i in range(start, _stop, stride):
        j = ':' + str(i + 1)
        # example: traj[':3']
        yield traj[j]


def find_lib(libname, unique=False):
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    lib_path_list = []
    key = "lib" + libname + "*"

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, key))
        for fname in fnamelist:
            if os.path.isfile(fname):
                lib_path_list.append(fname)

    if not lib_path_list:
        return None
    else:
        if unique:
            return set(lib_path_list)
        else:
            return lib_path_list


def read_orca_trj(fname):
    """return numpy 2D array
    """
    # http://stackoverflow.com/questions/14645789/
    #numpy-reading-file-with-filtering-lines-on-the-fly
    import numpy as np
    regexp = r'\s+\w+' + r'\s+([-.0-9]+)' * 3 + r'\s*\n'
    return np.fromregex(fname, regexp, dtype='f')


def read_gaussian_output(filename=None, top=None):
    """return a `pytraj.api.Trajectory` object

    Parameters
    ----------
    fname : str, filename
    top : {str, Topology}, optional, default None
        pytraj.Topology or a filename or None
        if None, use `antechamber` to generate mol2 file, need set $AMBERHOME env

    Requires
    --------
    cclib (``pip install cclib``)

    >>> import pytraj as pt
    >>> pt.tools.read_gaussian_output("gau.out", "mytest.pdb")
    """
    import pytraj as pt
    import cclib
    from pytraj.api import Trajectory
    from pytraj.utils.context import goto_temp_folder
    from pytraj._get_common_objects import _get_topology

    _top = _get_topology(None, top)
    gau = cclib.parser.Gaussian(filename)
    go = gau.parse()

    if _top is None:
        try:
            amberhome = os.environ['AMBERHOME']
        except KeyError:
            raise KeyError("must set AMBERHOME")

        fpath = os.path.abspath(filename)

        with goto_temp_folder():
            at = amberhome + "/bin/antechamber"
            out = "-i %s -fi gout -o tmp.mol2 -fo mol2 -at amber" % fpath
            cm = " ".join((at, out))
            os.system(cm)

            return Trajectory(xyz=go.atomcoords, top="tmp.mol2")
    else:
        return Trajectory(xyz=go.atomcoords, top=_top)


def read_to_array(fname):
    '''read text from file to numpy array'''
    import numpy as np
    with open(fname, 'r') as fh:
        arr0 = np.array([[x for x in line.split()] for line in fh.readlines()])
        return np.array(flatten(arr0), dtype='f8')


def merge_trajs(traj1, traj2, start_new_mol=True, n_frames=None):
    """

    Examples
    --------
       >>> # from two Trajectory or TrajectoryIterator
       >>> traj3 = merge_trajs(traj1, traj2)
       >>> assert traj3.n_frames == traj1.n_frames == traj2.n_frames
       >>> assert traj3.n_atoms == traj1.n_atoms + traj2.n_atoms
       >>> import numpy as np
       >>> assert np.any(traj3.xyz, np.vstack(tra1.xyz,  traj2.xyz)) == True

       >>> # from frame_iter for saving memory
       >>> traj3 = merge_trajs((traj1(0, 10, 2), traj1.top), 
                           (traj2(100, 110, 2), traj2.top), n_frames=6)

    Notes
    -----
    Code might be changed
    """
    from pytraj.compat import zip
    from pytraj import Trajectory
    import numpy as np

    if isinstance(traj1, (list, tuple)):
        n_frames_1 = n_frames
        top1 = traj1[1]
        _traj1 = traj1[0]
    else:
        n_frames_1 = traj1.n_frames
        top1 = traj1.top
        _traj1 = traj1

    if isinstance(traj2, (list, tuple)):
        n_frames_2 = n_frames
        top2 = traj2[1]  # example: (traj(0, 5), traj.top)
        _traj2 = traj2[0]
    else:
        n_frames_2 = traj2.n_frames
        top2 = traj2.top
        _traj2 = traj2

    if n_frames_1 != n_frames_2:
        raise ValueError("must have the same n_frames")

    traj = Trajectory()
    traj._allocate(n_frames_1, top1.n_atoms + top2.n_atoms)

    # merge Topology
    top = top1.copy()
    if start_new_mol:
        top.start_new_mol()
    top.join(top2)
    traj.top = top

    # update coords
    for f1, f2, frame in zip(_traj1, _traj2, traj):
        frame.xyz = np.vstack((f1.xyz, f2.xyz))

    return traj


def isel(traj, func, *args, **kwd):
    """iter-select frame based on func
    """
    for f in traj:
        if func(f, *args, **kwd):
            yield f
        else:
            pass


def filter(iterable, func):
    '''return a list
    '''
    return list(filter(func, iterable))


def as_2darray(traj):
    '''reshape traj.xyz to 2d array, shape=(n_frames, n_atoms * 3)

    Notes
    -----
    if ``traj`` is mutable, this method return a view of its coordinates.
    '''
    return traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)


def as_3darray(xyz):
    '''reshape xyz to 3d array, shape=(n_frames, n_atoms, 3)
    '''
    shape = xyz.shape
    if len(shape) != 2:
        raise ValueError('shape must be 2')
    new_shape = (shape[0], int(shape[1] / 3), 3)
    return xyz.reshape(new_shape)

def split_and_write_traj(self,
                         n_chunks=None,
                         root_name="trajx",
                         ext='nc', *args, **kwd):

    chunksize = self.n_frames // n_chunks
    for idx, traj in enumerate(self.iterchunk(chunksize=chunksize)):
        fname = ".".join((root_name, str(idx), ext))
        traj.save(fname, *args, **kwd)
