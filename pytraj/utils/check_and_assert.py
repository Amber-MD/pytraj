def file_exist(filename):
    import os
    return os.path.isfile(filename)

def is_generator(iter_obj):
    # use this method instead of `inspect` in python since this does not work with Cython
    # Reason: unknow
    if iter_obj.__class__.__name__ == 'generator':
        return True
    else:
        return False

def make_sure_exist(filename):
    if not file_exist(filename):
        txt = "can not find %s" % filename
        raise RuntimeError(txt)

def assert_almost_equal(arr0, arr1, decimals=3):
    '''numpy-like assert'''

    almost_equal = True
    SMALL = 10**(-decimals)

    for x, y in zip(arr0, arr1):
        if abs(x - y) > SMALL:
            almost_equal = False
    assert almost_equal == True

def _import_numpy():
    has_numpy = False
    try:
        __import__('numpy')
        has_numpy = True
        return (has_numpy, __import__('numpy'))
    except ImportError:
        has_numpy = False
        return (has_numpy, None)

def _import_h5py():
    has_h5py = False
    try:
        __import__('h5py')
        has_h5py = True
        return (has_h5py, __import__('h5py'))
    except ImportError:
        has_h5py = False
        return (has_h5py, None)

def _import(modname):
    """has_numpy, np = _import('numpy')"""
    has_module = False
    try:
        imported_mod = __import__(modname)
        has_module = True
        return (has_module, imported_mod)
    except ImportError:
        has_module = False
        return (has_module, None)

if __name__ == "__main__":
    import numpy as np
    assert_almost_equal([1., 2., 3.], [1., 2., 3.], decimals=3)
    assert_almost_equal([1., 2., 3.], [1., 2.,], decimals=3)
    #assert_almost_equal([1., 2., 4.], [1., 2., 3.], decimals=3)

    arr0 = np.array([1., 2., 3.])
    arr1 = np.array([1., 2., 3.])
    assert_almost_equal(arr0, arr1)
