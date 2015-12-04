from __future__ import absolute_import
import os
import numbers
from ..externals.six import string_types, zip


def eq(arr0, arr1):
    assert arr0 == arr1


def file_exist(filename):
    import os
    return os.path.isfile(filename)


def is_word_in_class_name(obj, word):
    """check if a `word` is in obj.__class__.__name__
    """
    return word in obj.__class__.__name__


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


def is_frame_iter(iter_obj):
    """check if is frame_iter
    """
    if iter_obj.__class__.__name__ == 'FrameIterator':
        return True
    return False


def is_int(num):
    """wrapper class to check if `num` is int
    isinstance(nu, (int, long)) does not work with numpy.int64, so we use numbers.Integral
    """
    return isinstance(num, numbers.Integral)


def is_number(num):
    return isinstance(num, numbers.Number)


def ensure_exist(filename):
    '''
    >>> ensure_exist('xfdasfda33fe')
    Traceback (most recent call last):
        ...
    RuntimeError: can not find xfdasfda33fe
    '''
    if not os.path.exists(filename):
        txt = "can not find %s" % filename
        raise RuntimeError(txt)


def ensure_not_none_or_string(obj):
    '''
    >>> ensure_not_none_or_string(None)
    Traceback (most recent call last):
        ...
    ValueError: <None> is a wrong input. Can not use `None` or string type
    >>> ensure_not_none_or_string('test')
    Traceback (most recent call last):
        ...
    ValueError: <test> is a wrong input. Can not use `None` or string type
    '''
    name = obj.__str__()
    msg = "<%s> is a wrong input. Can not use `None` or string type" % name
    if obj is None or isinstance(obj, string_types):
        raise ValueError(msg)


def assert_almost_equal(arr0, arr1, decimal=4):
    '''numpy-like assert,
    use default decimal=4 to match cpptraj's output
    >>> assert_almost_equal(0, 0)
    >>> assert_almost_equal([1, 2], [1.0000000003, 2.00000003])
    '''
    import math

    if is_number(arr0):
        arr0 = [arr0, ]
    if is_number(arr1):
        arr1 = [arr1, ]

    almost_equal = True
    SMALL = 10**(-decimal)

    if hasattr(arr0, 'flatten') and hasattr(arr1, 'flatten'):
        _arr0 = arr0.flatten()
        _arr1 = arr1.flatten()
    else:
        _arr0 = arr0
        _arr1 = arr1
    assert len(_arr0) == len(_arr1), 'two arrays must have the same length'

    for x, y in zip(_arr0, _arr1):
        if math.isnan(x) or math.isnan(y):
            raise ValueError('do not support NAN comparison')
        if abs(x - y) > SMALL:  # pragma: no cover
            almost_equal = False
    assert almost_equal is True


def _import(modname):
    """has_numpy, np = _import('numpy')
    >>> has_np, np = _import('numpy')
    """
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
    >>> has_np = has_("numpy")
    """
    return _import(lib)[0]


if __name__ == "__main__":
    import numpy as np
    assert_almost_equal([1., 2., 3.], [1., 2., 3.], decimals=3)
    assert_almost_equal([1., 2., 3.], [1., 2., ], decimals=3)

    arr0 = np.array([1., 2., 3.])
    arr1 = np.array([1., 2., 3.])
    assert_almost_equal(arr0, arr1)
