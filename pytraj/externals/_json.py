from __future__ import absolute_import
from .six import PY3


def to_json(obj, path):
    """
    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    import json
    if PY3:
        mode = 'w'
    else:
        mode = 'wb'

    with open(path, mode) as f:
        json.dump(obj, f)


def read_json(path):
    """
    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    dict : python dict
    """
    import json

    if PY3:
        mode = 'r'
    else:
        mode = 'rb'

    with open(path, mode) as fh:
        return json.load(fh)
