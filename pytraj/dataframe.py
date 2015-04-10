from pytraj.utils import _import
from pytraj import DataSetList

has_pandas, pd = _import("pandas")

if not has_pandas:
    print ("must have pandas installed to use DataFrame")

def to_dataframe(dslist):
    """return panda's DataFrame object
    Parameters
    ----------
    dslist : DataSetList object or list/tuple of DataSet
        or any object DataFrame can read
    """
    if isinstance(dslist, (list, tuple, DataSetList)):
        my_dict = dict((d0.legend, list(d0.to_ndarray())) for d0 in dslist)
    else:
        my_dict = dslist

    # what else I can do here?
    if has_pandas:
        return pd.DataFrame(my_dict)
    else:
        return None
