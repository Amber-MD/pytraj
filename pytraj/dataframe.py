from pytraj.utils import _import

has_pandas, pd = _import("pandas")

if not has_pandas:
    print ("must have pandas installed to use DataFrame")

def to_dataframe(dslist):
    """return panda's DataFrame object
    Parameters
    ----------
    dslist : DataSetList object
    """

    my_dict = {}

    for legend in dslist.get_legends():
        my_dict[legend] = dslist[legend][0][:]

    # what else I can do here?
    if has_pandas:
        return pd.DataFrame(my_dict)
    else:
        return None
