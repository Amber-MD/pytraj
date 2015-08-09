import os
import unittest
import numpy as np
from pytraj.base import *

ts = TrajectoryIterator()
ts.top = Topology("./data/Tc5b.top")
mdx = "./data/md1_prod.Tc5b.x"
ts.load(mdx)
print(ts.size)
farray = ts[:]
print(farray.size)


def tease_Traj(N=10, ts=ts):
    def forward():
        for i in range(N):
            print(ts[2])
            print(ts.total_read_frames)
            print(ts.current_frame)
            print(ts[2:])
            print(ts.total_read_frames)
            print(ts.current_frame)
            #assert ts[0:9:1][0].n_atoms == 304
            #assert ts[0:9][0].n_atoms == 304
            print(ts[0:9][0])
            print(ts.total_read_frames)
            print(ts.current_frame)
            ts1 = ts[0:9] + ts[0:9]
            print(ts.total_read_frames)
            print(ts.current_frame)
            print(ts1)

    def backward():
        for i in range(N):
            print(ts[-2])
            print(ts.total_read_frames)
            print(ts.current_frame)
            print(ts[-2][0])
            print(ts.total_read_frames)
            print(ts.current_frame)
            #assert ts[-1:-9:-1][0].n_atoms == 304
            print(ts[1:9:1][0])
            print(ts.total_read_frames)
            print(ts.current_frame)
            print(ts[-1:-9:-1][0])
            print(ts.total_read_frames)
            print(ts.current_frame)

    def mix():
        for i in range(N):
            print(ts[-2:-8:-1][0])
            ts[0:9:1][0] += ts[0:9:1][0]
            print(ts)

    # forward()
    # backward()
    mix()


tease_Traj(1, farray)
