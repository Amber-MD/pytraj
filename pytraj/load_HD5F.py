from pytraj.FrameArray import FrameArray
from pytraj.Frame import Frame

def load_hd5f(filename):
    farray = FrameArray()
    try:
        import h5py
        import numpy as np
    except ImportError:
        raise ImportError("require h5py, HD5F lib and numpy")

    fh = h5py.File(filename, 'r')
    crd = fh['coordinates']
    shape = crd.shape

    farray.resize(shape[0])
    for idx, arr in enumerate(crd):
        farray[idx] = Frame(shape[1])
        farray[idx].set_from_crd(arr.flatten().astype(np.float64))
    return farray
