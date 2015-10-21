from __future__ import print_function, absolute_import
import os
import warnings
from functools import wraps


# we duplicate code from .utils.check_and_assert here to avoid circular import

def _register_pmap(f):
    @wraps(f)
    def inner(*args, **kwd):
        return f(*args, **kwd)

    inner._is_parallelizable = True
    return inner

def _register_openmp(f):
    @wraps(f)
    def inner(*args, **kwd):
        return f(*args, **kwd)

    inner._openmp_capability = True
    return inner

def has_(lib):
    """check if having `lib` library

    Examples
    --------
    >>> has_("numpy")
    """
    return _import(lib)[0]


def makesureABC(classname):
    def inner(func):
        def _inner(self, *args, **kwd):
            if self.__class__.__name__ == classname:
                raise NotImplementedError("This is Abstract Base Class")
            else:
                return func(self, *args, **kwd)

        # update _inner doc for func
        _inner.__doc__ = func.__doc__
        _inner.__name__ = func.__name__
        return _inner

    return inner

def local_test(ext='edu'):
    import platform
    e = platform.node().split(".")[-1]
    if e != ext:
        do_test = False
    else:
        do_test = True

    def inner(func):
        def _no_test(*args, **kwd):
            if do_test:
                return func(*args, **kwd)
            else:
                txt = "skip. Only test on local"
                print(txt)
                return None

        return _no_test

    return inner

def deprecated(func):
    # from: https://wiki.python.org/moin/PythonDecoratorLibrary
    '''This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.'''

    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)

    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    return new_func
