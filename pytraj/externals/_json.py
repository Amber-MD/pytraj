from __future__ import absolute_import
import json
from .six import PY3

# adapted from pandas package
# see license in $PYTRAJHOME/license/externals/

def to_json(obj, path):
    """
    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
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

    if PY3:
        mode = 'r'
    else:
        mode = 'rb'

    with open(path, mode) as fh:
        return json.load(fh)
