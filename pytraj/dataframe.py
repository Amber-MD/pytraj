from __future__ import absolute_import
from .utils import _import, require
from .DataSetList import DataSetList


def to_dataframe(dslist):
    """return panda's DataFrame object
    Parameters
    ----------
    dslist : DataSetList object or list/tuple of DataSet
        or any object DataFrame can read
    """
    has_pandas, pd = _import_pandas()
    
    if not has_pandas:
        require("pandas")

    if isinstance(dslist, (list, tuple, DataSetList)):
        my_dict = dict((d0.legend, list(d0.to_ndarray())) for d0 in dslist)
    else:
        my_dict = dslist

    return pd.DataFrame(my_dict)
