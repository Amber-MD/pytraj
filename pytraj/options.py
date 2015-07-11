"""commonly used options
"""

# for TrajectoryIterator
# import pytraj as pt
# traj(**pt.options.iter_option)
iter_option = {
    'start': 0,
    'stop': None,
    'autoimage': True,
    'rmsfit': (0, '@CA')
}
