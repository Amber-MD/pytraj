from __future__ import absolute_import


def simple_plot(d0, *args, **kwd):
    from matplotlib import pyplot as plt
    fig = plt.pyplot.plot(range(d0.size), d0[:], *args, **kwd)
    return fig
