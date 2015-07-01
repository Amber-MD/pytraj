from __future__ import absolute_import
import sys
from itertools import islice, groupby
import functools
from collections import OrderedDict
from pytraj.externals.six import iteritems

# string_types, PY2, PY3 is copied from six.py
# see license in $PYTRAJHOME/license/externals/
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str
else:
    string_types = basestring

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


def split_range(n_chunks, start, stop):
    '''
    >>> from pytraj.misc import split_range
    >>> split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    list_of_tuple = []
    chunksize = (stop - start) // n_chunks
    for i in range(n_chunks):
        if i < n_chunks - 1:
            _stop = (i + 1) * chunksize
        else:
            _stop = stop
        list_of_tuple.append((start + i * chunksize, _stop))
    return list_of_tuple


@_dispatch_value
def split(data, n_chunks_or_array):
    """split `self.data` to n_chunks

    Notes : require numpy (same as `array_split`)
    """
    return np.array_split(data, n_chunks_or_array)


def chunk_average(self, n_chunk, restype='same'):
    import numpy as np
    data = np.array(list(map(np.mean, split(self, n_chunk))))
    if restype == 'same':
        new_array = self.shallow_copy()
        new_array.values = data
        return new_array
    else:
        return data


@_not_yet_tested
def moving_average(data, n):
    # http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    """
    Note
    ----
    not assert yet
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
    # copied from pandas PR
    # https://github.com/ghl3/pandas/blob/groupby-pipe/pandas/tools/util.py
    # see license in pytraj/license/
    """
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
    # copied from pandas (added pytraj's example)
    # see license in pytraj/license/
    """Compose 2 or more callables

    Examples
    --------
    >>> import pytraj as pt
    >>> func = pt.tools.compose(pt.calc_radgyr, pt.iterload)
    >>> func("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    """
    assert len(funcs) > 1, 'At least 2 callables must be passed to compose'
    return reduce(_compose2, funcs)

def grep(self, key):
    """
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
    # http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
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


@_not_yet_tested
def n_grams(a, n):
    """
    """
    # http://sahandsaba.com/thirty-python-language-features-and-tricks-you-may-not-know.html
    z = (islice(a, i, None) for i in range(n))
    return zip(*z)


def dict_to_ndarray(dict_of_array):
    """
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
    """
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
    return np.vstack((np.array(f.xyz) for f in iterables))


def merge_frames(iterables):
    """merge_coordinates from frames
    """
    from pytraj import Frame
    xyz = np.vstack((np.array(f.xyz) for f in iterables))
    return Frame().append_xyz(xyz)


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

    tmp = sum((arr1 - arr2)**2)
    return sqrt(tmp / arr1.shape[0])


def rmsd(a1, a2, flatten=True):
    """
    rmsd for two array with the same shape

    Parameters
    ----------
    a1, a2: np.ndarray
    """
    if a1.shape != a2.shape and not flatten:
        raise ValueError("must have the same shape")
    return rmsd_1darray(a1.flatten(), a2.flatten())
