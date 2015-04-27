"""simple plot"""
from __future__ import absolute_import
from pytraj.utils.check_and_assert import require
from .plot_matrix import plot_matrix

def show_config():
    """show good ipython config"""
    txt = """
    %matplotlib inline # inline for matplotlibe
    %config InlineBackend.figure_format = 'retina'  # high resolution
    import matplotlib
    matplotlib.rcParams['savefig.dpi'] = 2 * matplotlib.rcParams['savefig.dpi'] # larger image
    """
    print (txt)
