# that stops travis
import os
import pytraj as pt
from numpy import max
from memory_profiler import profile, memory_usage

from glob import glob

traj0 = pt.load_sample_data('tz2')
# make fake large trajectory from 10 frames
for _ in range(9):
    traj0.load(traj0.filelist)

tlist = [traj0, ]
if os.path.exists("./data/nogit/remd/myparm.parm7"):
    tlist.append(pt.iterload(
        glob("./data/nogit/remd/remd.x.*")[:10],
        "./data/nogit/remd/myparm.parm7"))

for traj in tlist[:1]:
    print(traj)
    print(traj._estimated_GB)

    @profile
    def test_pairwise_rmsd():
        pt.pairwise_rmsd(traj(stop=1000), mask='@CA')

    @profile
    def test_write():
        pt.write_traj('test.nc',
                      traj,
                      frame_indices=range(10000),
                      overwrite=True)
        print(pt.iterload('./test.nc', traj.top))

    @profile
    def test():
        # OK: no memory problem
        for frame in traj:
            pass

    @profile
    def test_closest():
        print('test_closest')
        for frame in pt.closest(
                traj(0, 100),
                n_solvents=10,
                restype='iterator')[0]:
            pass

    @profile
    def test_strip(chunksize=100):
        print("test_strip")
        for chunk in traj.iterchunk(chunksize=100):
            chunk.strip(':WAT')

    @profile
    def test_iterchunk(chunksize=100):
        # OK
        for chunk in traj.iterchunk(chunksize=chunksize):
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
    def test_autoimage_iterchunk_0(chunksize=100):
        # OK
        for chunk in traj.iterchunk(chunksize=chunksize):
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
            pt.translate(frame, '', top=traj.top)

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

    @profile
    def test_iter_frame_indices():
        print('test_iter_frame_indices')
        for idx, f in enumerate(traj(frame_indices=range(0, traj.n_frames,
                                                         3))):
            pass
        print(idx)

    func_list = [
        test_iterchunk,  # test_pairwise_rmsd,
        # test_write,
        # test_strip,
        # test_closest,
        # test_iter_frame_indices,
        # test_center,
        # test_rmsd,
        # test_simple_frame_iter,
        # test_frame_iter_with_mask,
        # test_autoimage,
        # test_autoimage_regular,
        # test_translate_regular,
        # test_iterchunk,
        # test_autoimage_iterchunk_0,
    ]

    estimated_GB = traj._estimated_GB

    for func in func_list:
        mem = max(memory_usage(func))
        print("%s : %s" % (func.__name__, mem))
