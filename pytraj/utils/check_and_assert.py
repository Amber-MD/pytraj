from __future__ import absolute_import
import os
import numbers
from ..externals.six import string_types
from functools import wraps

# don't `import pytraj.externals.six` here: got import error
try:
    zip = izip
except:
    pass


def eq(arr0, arr1):
    assert arr0 == arr1


def eq_coords(fa0, fa1):
    # use this method for very large trajs
    # use `assert_almost_equal` for `xyz` is to slow since need to copy to xyz
    # first
    count = 0
    import numpy as np
    for f0, f1 in zip(fa0, fa1):
        count += 1
        assert np.any(f0.xyz == f1.xyz) == True
    assert count == fa0.n_frames == fa1.n_frames


def a_isinstance(obj, class_type):
    assert isinstance(obj, class_type) == True


def file_exist(filename):
    import os
    return os.path.isfile(filename)


def get_amber_saved_test_dir(suffix):
    """return full dir of amber test file or None
    Used for assert in testing

    Parameter
    --------
    suffix : str
    """
    import os
    try:
        amberhome = os.environ['AMBERHOME']
        return amberhome + "/AmberTools/test/cpptraj/" + suffix
    except:
        return None


def is_linux():
    import sys
    return 'linux' in sys.platform


def is_word_in_class_name(obj, word):
    """check if a `word` is in obj.__class__.__name__
    """
    return word in obj.__class__.__name__


def is_pytraj_trajectory(obj):
    return is_word_in_class_name(obj, 'Trajectory')


def is_range(obj):
    return is_word_in_class_name(obj, 'range')


def is_array(obj):
    """check if a `word` is in obj.__class__.__name__
    """
    if is_word_in_class_name(obj, 'array'):
        return True
    else:
        return False


def are_instance(obj_list, cls):
    """check if all elements have the same class `cls`"""
    for element in obj_list:
        if not isinstance(element, cls):
            return False
    return True


def is_generator(iter_obj):
    # use this method instead of `inspect` in python since this does not work with Cython
    # Reason: (I don't know)
    if iter_obj.__class__.__name__ == 'generator':
        return True
    else:
        return False


def is_mdanalysis(obj):
    return is_word_in_class_name(obj, 'Universe')


def is_mdtraj(obj):
    """check if traj is mdtraj object"""
    return True if 'mdtraj' in obj.__str__() else False


def is_mdanalysis(obj):
    return is_word_in_class_name(obj, 'Universe')


def is_frame_iter(iter_obj):
    """check if is frame_iter

    See Also
    --------
    Trajectory.frame_iter
    Trajin.frame_iter
    """
    if iter_obj.__class__.__name__ == 'generator' and 'frame_iter' in iter_obj.__name__:
        return True
    if iter_obj.__class__.__name__ == 'FrameIter':
        return True
    return False


def is_chunk_iter(iter_obj):
    """check if is frame_iter

    See Also
    --------
    Trajectory.frame_iter
    Trajin.frame_iter
    """
    try:
        name = iter_obj.__name__
    except AttributeError:
        return False

    if iter_obj.__class__.__name__ == 'generator' and (
            'chunk_iter' in name or 'iterchunk' in name):
        return True
    else:
        return False


def is_int(num):
    """wrapper class to check if `num` is int
    isinstance(nu, (int, long)) does not work with numpy.int64, so we use numbers.Integral
    """
    return isinstance(num, numbers.Integral)


def is_number(num):
    return isinstance(num, numbers.Number)


def ensure_exist(filename):
    if not os.path.exists(filename):
        txt = "can not find %s" % filename
        raise RuntimeError(txt)


def ensure_not_none_or_string(obj):
    name = obj.__str__()
    msg = "<%s> is a wrong input. Can not use `None` or string type" % name
    if obj is None or isinstance(obj, string_types):
        raise ValueError(msg)


def assert_almost_equal(arr0, arr1, decimal=3):
    '''numpy-like assert'''

    if is_number(arr0):
        arr0 = [arr0, ]
    if is_number(arr1):
        arr1 = [arr1, ]

    almost_equal = True
    SMALL = 10 ** (-decimal)

    if hasattr(arr0, 'flatten') and hasattr(arr1, 'flatten'):
        _arr0 = arr0.flatten()
        _arr1 = arr1.flatten()
    else:
        _arr0 = arr0
        _arr1 = arr1
    assert len(_arr0) == len(_arr1), 'two arrays must have the same length'

    for x, y in zip(_arr0, _arr1):
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


def _import_pandas():
    has_pd = False
    try:
        pd = __import__('pandas')
        has_pd = True
        # set print options
        pd.options.display.max_rows = 20
        return (has_pd, pd)
    except ImportError:
        has_pd = False
        return (has_pd, None)


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


def has_(lib):
    """check if having `lib` library
    Example:
    >>> has_("numpy")
    """
    return _import(lib)[0]


def require(libname):
    has_lib, lib = _import(libname)
    if not has_lib:
        txt = "require %s lib" % libname
        raise ImportError(txt)


if __name__ == "__main__":
    import numpy as np
    assert_almost_equal([1., 2., 3.], [1., 2., 3.], decimals=3)
    assert_almost_equal([1., 2., 3.], [1., 2., ], decimals=3)
    #assert_almost_equal([1., 2., 4.], [1., 2., 3.], decimals=3)

    arr0 = np.array([1., 2., 3.])
    arr1 = np.array([1., 2., 3.])
    assert_almost_equal(arr0, arr1)
