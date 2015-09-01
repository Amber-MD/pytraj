def joinplot(data, x='', y='', show=True, **kwd):
    """require seaborn, pandas

    See also
    --------
    seaborn.jointplot
    """
    import seaborn as sb
    from matplotlib import pyplot as plt

    if not x:
        x = data[0].key
    if not y:
        y = data[1].key

    ax = sb.jointplot(x=x, y=y, data=data.to_dataframe(), **kwd)
    if show:
        plt.show()
    return ax
