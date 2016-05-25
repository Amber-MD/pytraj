#!/usr/bin/env python
import numpy as np
import pytraj as pt
from memory_profiler import profile
from pytraj.testing import get_fn

fn, tn = get_fn('tz2')

traj = pt.iterload([fn, ]*100, tn)
print(traj)
ref = traj[:1]

@profile
def do_it(traj=traj,n_times=1):

    # memory leaking
    traj.superpose(ref=ref, mask='@CA')

    # memory not leaking
    # traj.autoimage().center()

    for _ in range(n_times):
        for frame  in traj: pass

    print(traj._cdslist)

if __name__ == '__main__':
    # expectation: memory usage should be identical/similiar if varying n_times

    import sys
    n_times = int(sys.argv[1])

    do_it(n_times=n_times)
