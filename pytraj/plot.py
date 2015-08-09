from __future__ import absolute_import
from .compat import string_types
from .utils import is_int


def plot(x='', y='', data=None, *args, **kwd):
    from matplotlib import pyplot as plt
    try:
        import seaborn as sb
        sb.set()
    except ImportError:
        pass
    ax = plt.subplot(111)
    if x == '':
        xd = data
    else:
        xd = data[x]

    if y != '':
        yd = data[y]
        ax.plot(xd, yd, *args, **kwd)
    else:
        ax.plot(xd, *args, **kwd)
    return ax


import sys
sys.modules[__name__] = plot
