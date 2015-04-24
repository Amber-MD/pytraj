from __future__ import absolute_import
from .utils import _import

def simple_plot(d0, *args, **kwd):
    has_plot, plt = _import('matplotlib.pyplot')
    if not has_plot:
        raise RuntimeError("require matplotlib installed")
    fig = plt.pyplot.plot(range(d0.size), d0[:], *args, **kwd)
    return fig
