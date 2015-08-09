'''pytraj note: changed (a bit) from pandas.io.pickle module
'''
# pandas is distributed under a 3-clause ("Simplified" or "New") BSD
# license.

# Note: full license is in $PYTRAJHOME/license/externals/pandas.txt

from __future__ import absolute_import
from .six.moves import cPickle as pkl
from .six import PY3


def to_pickle(obj, path):
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    with open(path, 'wb') as f:
        pkl.dump(obj, f, protocol=pkl.HIGHEST_PROTOCOL)


def read_pickle(path):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Warning: Loading pickled data received from untrusted sources can be
    unsafe. See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """

    def try_read(path, encoding=None):
        # try with cPickle
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param

        # cpickle
        # GH 6899
        try:
            with open(path, 'rb') as fh:
                return pkl.load(fh)
        except (Exception) as e:

            # reg/patched pickle
            try:
                with open(path, 'rb') as fh:
                    return pc.load(fh, encoding=encoding, compat=False)

            # compat pickle
            except:
                raise NotImplementedError("don't know how to read")
                # with open(path, 'rb') as fh:
                #    return pc.load(fh, encoding=encoding, compat=True)

    try:
        return try_read(path)
    except:
        if PY3:
            return try_read(path, encoding='latin1')
        raise
