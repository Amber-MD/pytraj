#!/usr/bin/env python

import pytraj as pt

traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo')

def test(traj=traj):
    return pt.molsurf(traj, '@P')

def test_2(traj=traj):
    data = []
    for idx, c in enumerate(pt.iterchunk(traj, 100)):
        c.strip_atoms('!@P')
        print(idx, c)
        data.append(pt.molsurf(c))
    return data

data = test()
#data = test_2()
print(data)
# result:
