"""this file has commonly used actions such as rmsd calculation, 
randomizeions, strip atoms, ..."""

from __future__ import print_function, absolute_import
import os
from glob import glob
from pytraj.cpp_options import set_world_silent
from pytraj.compat import set

# external
from pytraj.externals.six import string_types

try:
    from pytraj.externals.magic import from_file as file_type_info
except ImportError:
    file_type_info = None


def info(obj=None):  # pragma: no cover
    """get `help` for obj
    Useful for Actions and Analyses

    Since we use `set_worl_silent` to turn-off cpptraj' stdout, we need 
    to turn on to use cpptraj's help methods
    """
    from pytraj import adict, analdict
    adict_keys = adict.keys()
    anal_keys = analdict.keys()

    if obj is None:
        print("action's keys", adict_keys)
        print("analysis' keys", anal_keys)
    else:
        if isinstance(obj, string_types):
            if obj in adict.keys():
                # make Action object
                _obj = adict[obj]
            elif obj in analdict.keys():
                # make Analysis object
                _obj = analdict[obj]
            else:
                raise ValueError("keyword must be an Action or Analysis")
        else:
            # assume `obj` hasattr `help`
            _obj = obj

        if hasattr(_obj, 'help'):
            set_world_silent(False)
            _obj.help()
            set_world_silent(True)
        elif hasattr(_obj, 'info'):
            set_world_silent(False)
            _obj.info()
            set_world_silent(True)
        elif 'calc_' in _obj.__name__:
            key = _obj.__name__.split("_")[-1]
            set_world_silent(False)
            adict[key].help()
            set_world_silent(True)
        elif hasattr(_obj, '__doc__'):
            print(_obj.__doc__)
        else:
            raise ValueError("object does not have `help` method")


def show_code(func, get_txt=False):  # pragma: no cover
    """show code of func or module"""
    import inspect
    txt = inspect.getsource(func)
    if not get_txt:
        print(txt)
    else:
        return txt


def get_atts(obj):  # pragma: no cover
    """get methods and atts from obj but excluding special methods __"""
    atts_dict = dir(obj)
    return [a for a in atts_dict if not a.startswith("__")]


def find_libcpptraj(**kwd):  # pragma: no cover
    '''
    '''
    return find_library('cpptraj', **kwd)


def find_library(libname, unique=False):  # pragma: no cover
    """return a list of all library files"""
    paths = os.environ.get('LD_LIBRARY_PATH', '').split(':')
    lib_path_list = []
    key = "lib" + libname + "*"

    for path in paths:
        path = path.strip()
        fnamelist = glob(os.path.join(path, key))
        for fname in fnamelist:
            if os.path.isfile(fname):
                lib_path_list.append(fname)

    if not lib_path_list:
        return None
    else:
        if unique:
            return set(lib_path_list)
        else:
            return lib_path_list