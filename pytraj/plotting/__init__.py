"""simple plot"""
from __future__ import absolute_import
from .plot_matrix import plot_matrix
from .wrap_seaborn import joinplot

try:
    from matplotlib.pyplot import show, plot
except ImportError:
    show = None
    plot = None


_pylab_config = """
%matplotlib inline # inline for matplotlibe
%config InlineBackend.figure_format = 'retina'  # high resolution
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi'] # larger image
"""

def show_config():
    """show good ipython config"""
    return _pylab_config
