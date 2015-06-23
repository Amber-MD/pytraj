from pytraj.utils import _import_numpy

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
