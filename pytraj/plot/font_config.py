'''config font, color, ... for plotting
'''
try:
    from matplotlib.pyplot import rc

    font = {'family': 'serif', 'size': '14'}
    rc('font', **font)
except ImportError:
    font = None
    rc = None
