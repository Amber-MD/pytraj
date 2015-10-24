"""simple plot"""
from __future__ import absolute_import
from .plot_matrix import plot_matrix
from .wrap_seaborn import joinplot
from . import symbols

try:
    from matplotlib.pyplot import rc
    from matplotlib.pyplot import show, plot
    from matplotlib import pyplot as plt
    font = {'family': 'serif', 'size': '15'}
    rc('font', **font)
except ImportError:
    show = None
    plot = None
    plt = None
    font = None
    rc = None

_pylab_config = """
%matplotlib inline # inline for matplotlibe
%config InlineBackend.figure_format = 'retina'  # high resolution
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi'] # larger image
"""


def show_config():
    """show good ipython config"""
    return _pylab_config


def polar(data, *args, **kwd):
    ax = plt.subplot(111, polar=True)
    ax.plot(data, *args, **kwd)
    return ax


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
