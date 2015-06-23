from pytraj.utils import _import_numpy
from pytraj.compat import string_types

try:
    # PY3
    from functools import reduce
except ImportError:
    # 
    pass

# this module gathers commonly used functions
# from toolz, stackoverflow, and from myself

_, np = _import_numpy()

def _dispatch_value(func):
    def inner(data, *args, **kwd):
        if hasattr(data, 'values'):
            _data = data.values
        else:
            _data = data
        return func(_data, *args, **kwd)
    inner.__doc__ = func.__doc__
    return inner

def _no_tested(func):
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

@_dispatch_value
def chunk_average(self, n_chunk):
    import numpy as np
    return np.array(list(map(np.mean, split(self, n_chunk))))

@_dispatch_value
@_no_tested
def moving_average(data, n):
    # http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    """
    Note
    ----
    not assert yet
    """
    window = np.ones(int(n))/float(n)
    return np.convolve(data, window, 'same')

@_dispatch_value
def pipe(self, *funcs):
    """apply a series of functions to self's data
    """
    if hasattr(self, 'values'):
        values = self.values
    else:
        values = self

    for func in funcs:
        values = func(values)
    return values

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
        #if isinstance(el, (list, tuple)):
        if hasattr(el, "__iter__") and not isinstance(el, string_types):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result
