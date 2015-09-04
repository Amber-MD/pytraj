# no unittest (and no travis test)
from __future__ import print_function
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.utils import Timer
from pytraj._shared_methods import _frame_iter
from pytraj.compat import zip


def test_0():
    traj = mdio.iterload("./data/tz2.ortho.nc", "./data/tz2.ortho.parm7")
    traj_saved = mdio.iterload("./data/tz2.autoimage_with_rmsfit.nc", traj.top)
    fa = traj[:]
    ref0 = fa[0].copy()
    ref1 = fa[0].copy()
    mask = '@H='

    # test autoimage
    for chunk in traj.chunk_iter(chunksize=2, autoimage=True):
        pass
    fa0 = fa.copy()
    fa0.autoimage()
    aa_eq(chunk.xyz, fa0[-2:].xyz)

    # test rmsfit
    for chunk in traj.chunk_iter(chunksize=2, rmsfit=(ref0, mask)):
        pass
    fa0 = fa.copy()
    fa0.rmsfit(ref0, mask)
    aa_eq(chunk.xyz, fa0[-2:].xyz)

    # test rmsfit and autoimage
    for chunk in traj.chunk_iter(chunksize=2,
                                 rmsfit=(ref0, mask),
                                 autoimage=True):
        pass

    fa0 = fa.copy()
    fa0.autoimage()
    fa0.rmsfit(ref0, mask)
    aa_eq(chunk.xyz, fa0[-2:].xyz)
    print(traj[-2:].xyz[0, 0], chunk.xyz[0, 0], fa0[-2:].xyz[0, 0])

    # assert to cpptraj: need to set mass
    fa1 = fa.copy()
    for frame in fa1:
        frame.set_frame_mass(fa1.top)
    fa1.autoimage()
    fa1.rmsfit(ref0, mask)
    print(traj_saved[-2:].xyz[0, 0], fa1[-2:].xyz[0, 0])
    # FIXME: assert failed
    #aa_eq(fa1.xyz, traj_saved[-2:].xyz)
    for f0, f1 in zip(traj_saved, fa1):
        print(f0.rmsd_nofit(f1), f0.rmsd(f1))

    fa_saved = traj_saved['!:WAT']
    fa1_nowat = fa1['!:WAT']
    # FIXME: assert failed
    #aa_eq(fa_saved.xyz, fa1_nowat.xyz)


if __name__ == "__main__":
    test_0()
