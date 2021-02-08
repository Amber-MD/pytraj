import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq


def test():
    traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
    itertraj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

    for idx, (f0, f1) in enumerate(zip(traj, itertraj)):
        aa_eq(f0.xyz, f1.xyz)
    assert idx == traj.n_frames - 1
