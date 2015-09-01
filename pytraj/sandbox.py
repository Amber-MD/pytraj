'''sandbox for Julialang and other stuff.
Codes in this module might not be tested and they might be not correct, at all.
'''


def take(traj, indices):
    return traj[indices]


def itake(traj, indices):
    return traj.iterframe(frame_indices=indices)


def get_top(traj):
    return traj.top


def set_top(traj, top):
    traj.top = top
