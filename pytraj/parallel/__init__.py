from .map import map
from pytraj.tools import concat_dict


def gather(name='data', clients=None, restype='ndarray'):
    '''gather data from different clients

    Parameters
    ----------
    name : name of the output holder
        for example: data = pytraj.calc_radgyr(traj) --> name = 'data'
    clients : IPython.parallel.Client objects
        number of clients == n_cores you use
    restype : str, {'ndarray', 'dataset'}, default 'ndarray' 
        if 'ndarray': hstack data by numpy.vstack
        if 'dataset': 'data' should be a list of dict, then will be converted
        to `pytraj.datasetlist.DatasetList` object

    Examples
    --------
    (fill me)
    '''
    if restype == 'ndarray':
        import numpy as np
        return np.hstack((x[name] for x in clients))
    elif restype == 'dataset':
        # it's user's responsibility to return a list of dicts
        from pytraj import datasetlist
        iter_of_dslist = (
            datasetlist._from_full_dict(x[name]) for x in clients)
        return datasetlist.vstack(iter_of_dslist)
    elif restype == 'dict':
        return concat_dict((x[name] for x in clients))
    else:
        raise ValueError("must be ndarray | dataset | dict")
