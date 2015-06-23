from pytraj.utils import _import_numpy

# this module gathers commonly used functions
# from toolz, stackoverflow, and from myself

_, np = _import_numpy()

def chunk_average(data):
    pass

def slicing_sum(data, windows):
    pass

def slicing_mean(data, windows):
    N = windows
    return np.convolve(data, np.ones((N,))/N)[(N-1):]

def slicing_max(data, windows):
    pass

def slicing_min(data, windows):
    pass

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
