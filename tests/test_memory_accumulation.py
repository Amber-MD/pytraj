#import unittest # turn off to avoid lengthy non-output
# that stops travis
import os
import pytraj as pt
from numpy import max
from memory_profiler import profile, memory_usage

from glob import glob

traj0 = pt.load_sample_data('tz2')
# make fake large trajectory from 10 frames
for _ in range(10):
    traj0.load(traj0.filelist)

tlist = [traj0, ]
if os.path.exists("./data/nogit/remd/myparm.parm7"):
    tlist.append(pt.iterload(glob("./data/nogit/remd/remd.x.*")[:10],
                             "./data/nogit/remd/myparm.parm7"))

for traj in tlist:
    print(traj)
    print(traj._estimated_MB)

    @profile
    def test():
        # OK: no memory problem
        for frame in traj:
            pass

    @profile
    def test_chunk_iter(chunksize=100):
        # OK
        for chunk in traj.chunk_iter(chunksize=chunksize):
            pass

    @profile
    def test_simple_frame_iter():
        # OK: no memory problem
        for frame in traj():
            pass

    @profile
    def test_frame_iter_with_mask():
        # OK: no memory problem
        for frame in traj(mask='@CA'):
            pass

    @profile
    def test_autoimage():
        # OK
        for frame in traj(autoimage=True):
            pass

    @profile
    def test_autoimage_chunk_iter_0(chunksize=100):
        # OK
        for chunk in traj.chunk_iter(chunksize=chunksize):
            chunk.autoimage()

    @profile
    def test_autoimage_regular():
        # OK
        for frame in traj(stop=400):
            pt.autoimage(frame, top=traj.top)

    @profile
    def test_translate_regular():
        # OK
        for frame in traj(stop=400):
            pt.common_actions.translate(frame, '', top=traj.top)

    @profile
    def test_rmsd():
        # OK
        pt.rmsd(traj, ref=0)

    @profile
    def test_center():
        # OK
        for f in traj:
            pt.center(f, top=traj.top)

    @profile
    def test_molsurf():
        # OK
        pt.molsurf(traj(stop=-3, stride=2, autoimage=True, rmsfit=0))

    func_list = [
        test_center,
        test_rmsd,
        test_simple_frame_iter,
        test_frame_iter_with_mask,
        test_autoimage,
        test_autoimage_regular,
        test_translate_regular,
        test_chunk_iter,
        test_autoimage_chunk_iter_0,
    ]

    estimated_MB = traj._estimated_MB

    for func in func_list:
        func()
