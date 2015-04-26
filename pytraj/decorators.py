from __future__ import print_function, absolute_import
from .utils import _import
import os

# we duplicate code from .utils.check_and_assert here to avoid circular import
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

def for_testing(func):
    def inner(*args, **kwd):
        print ("this %s method is for tesing purpose" % func.__name__ )
        return func(*args, **kwd)
    return inner

def iter_warning(func):
    def inner(*args, **kwd):
        if args[0].size <= 0:
            raise ValueError("empty object, cannot do iteration")
        return func(*args, **kwd)
    return inner

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

def name_will_be_changed(msg):
    txt = "this method name will be changed"
    if not msg == "":
        txt +=  "to %s" % msg
    def inner(func):
        def _inner(self, *args, **kwd):
            print (txt)
            return func(self, *args, **kwd)
        return _inner
    return inner

def no_test(func):
    def _no_test(*args, **kwd):
        pass
    return _no_test

def test_if_having(lib):
    def inner(func):
        def _no_test(*args, **kwd):
            if has_(lib):
                return func(*args, **kwd)
            else:
                txt = "Does not have %s. Skip test" % lib
                print(txt)
                return None
        return _no_test
    return inner

def test_if_path_exists(mydir):
    def inner(func):
        def _no_test(*args, **kwd):
            if os.path.exists(mydir):
                return func(*args, **kwd)
            else:
                txt = "%s does not exist. Skip test" % mydir
                print(txt)
                return None
        return _no_test
    return inner

def require_having(mylib):
    def inner(func):
        def more_inner(*args, **kwd):
            has_lib, lib = _import(mylib)
            if not has_lib:
                msg = "need %s" % mylib
                raise ImportError(msg)
            else:
                return func(*args, **kwd)
        return more_inner
    return inner

def not_yet_supported(func):
    def inner(*args, **kwd):
        print ("%s not_yet_supported" % func.__name__)
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
