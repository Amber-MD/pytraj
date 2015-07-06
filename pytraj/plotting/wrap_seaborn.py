
def joinplot(data, x='', y='', **kwd):
    """require seaborn, pandas

    See also
    --------
    seaborn.jointplot
    """
    import seaborn as sb

    if not x:
        x = data[0].key
    if not y:
        y = data[1].key

    return sb.jointplot(x=x, y=y, data=data.to_dataframe(), **kwd) 
