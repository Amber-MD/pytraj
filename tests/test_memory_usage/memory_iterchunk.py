# that stops travis
import pytraj as pt
from utils import fn
from memory_profiler import profile


traj = pt.load_sample_data('tz2')
# make fake large trajectory from 10 frames
for _ in range(9):
    traj.load(traj.filelist)


@profile
def test_iterchunk(chunksize=100):
    for chunk in traj.iterchunk(chunksize=chunksize):
        pass


test_iterchunk(chunksize=100)
