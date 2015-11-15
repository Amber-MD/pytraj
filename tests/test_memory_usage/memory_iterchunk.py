# that stops travis
import os
import pytraj as pt
from numpy import max
from memory_profiler import profile, memory_usage

from glob import glob

traj = pt.load_sample_data('tz2')
# make fake large trajectory from 10 frames
for _ in range(9):
    traj.load(traj.filelist)


@profile
def test_iterchunk(chunksize=100):
    for chunk in traj.iterchunk(chunksize=chunksize):
        pass


test_iterchunk(chunksize=100)
