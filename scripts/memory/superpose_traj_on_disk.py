#!/usr/bin/env python
import pytraj as pt
from memory_profiler import profile
from pytraj.testing import get_fn

fn, tn = get_fn('tz2')

traj = pt.iterload([fn, ]*100, tn)

@profile
def do_it(traj=traj):
    ref = traj[:1]

    traj.superpose(ref=ref, mask='@CA')

    for _ in traj: pass
    print(traj._cdslist['__myrmsd'].values.shape)

    for _ in traj: pass
    print(traj._cdslist['__myrmsd'].values.shape)

    for _ in traj: pass
    print(traj._cdslist['__myrmsd'].values.shape)

    for _ in traj: pass
    print(traj._cdslist['__myrmsd'].values.shape)

    for _ in traj: pass
    print(traj._cdslist['__myrmsd'].values.shape)

do_it()
