import json

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
    with open(path, 'wb') as f:
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

    def try_read(path):
        try:
            with open(path, 'rb') as fh:
                return json.load(fh)
        except (Exception) as e:
            raise
    return try_read(path)
