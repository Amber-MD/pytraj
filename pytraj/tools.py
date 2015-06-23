from pytraj.utils import _import_numpy

_, np = _import_numpy()

def chunk_average(data):
    pass

def slicing_sum(data, windows):
    pass

def slicing_mean(data, windows):
    N = windows
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def slicing_max(data, windows):
    pass

def slicing_min(data, windows):
    pass
