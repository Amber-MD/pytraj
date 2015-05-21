from __future__ import print_function
#import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.testing import make_fake_traj
import pytraj.common_actions as pyca
from memory_profiler import memory_usage
from pytraj.compat import range
import time

MAXMEM = 1275 # MB
MINMEM = 1370
MAX_FRAME_MEM = 15. # MB
#XYZ_MEM = 
n_frames = 800
n_atoms = 50000
print ("MAXMEM=%s, MINMEM=%s, MAX_FRAME_MEM=%s" % (MAXMEM, MINMEM, MAX_FRAME_MEM))

traj = make_fake_traj(n_frames, n_atoms)

def single_traj():
    traj

def traj_xyz():
    traj.xyz

def single_frame():
    frame = traj[0].copy()

def make_random_traj():
    traj = make_fake_traj(n_frames, n_atoms)

def slicing():
    # make sure that slicing with copy=False does not  introduce much memory
    traj[:20]
    traj[:]
    traj.join([traj[:], traj[:], traj[:]], copy=False)

def iterating():
    for frame in traj:
        pass

def clustering():
    pyca.do_clustering(traj[:100], "kmeans clusters 2")

def make_copy_2_frames():
    f0 = traj[0].copy()
    #f1 = traj[1].copy()

def perform_action_radgyr():
    pyca.calc_radgyr(traj)

def perform_action_radgyr_dset_view():
    ds = pyca.calc_radgyr(traj, dtype='dataset')
    d = ds[:]

def perform_action_radgyr_dset_copy():
    ds = pyca.calc_radgyr(traj, dtype='dataset')
    d = ds.copy()

def inline_math_add(traj):
    traj += 1.

def apply_func(traj):
    traj.apply(lambda x : x * 2)

if __name__ == "__main__":
    import numpy as np
    from numpy import max

    m = np.max(memory_usage(single_traj))
    print ("single_traj", m)
    real_m_traj = m

    m = max(memory_usage((inline_math_add, (traj,))))
    print ('inline_math_add', m - real_m_traj)
    assert (m - real_m_traj) < 2.

    m = max(memory_usage((apply_func, (traj,))))
    print ('apply_func', m - real_m_traj)

    m = memory_usage(single_frame)[-1]
    print ("single_frame", m - real_m_traj)
    real_m_frame = m
    assert (m - MAX_FRAME_MEM  < MAXMEM)

    m = max(memory_usage(traj_xyz))
    print ("traj_xyz", m)
    assert (m - MAXMEM) < MAXMEM
    m_top = real_m_traj - (m - real_m_traj)
    print ("m_top", m_top)

    m = memory_usage(slicing)[-1]
    print (m)
    assert m - real_m_traj < real_m_frame

    m = max(memory_usage(iterating))
    print ('iterating', m - real_m_traj)
    assert m - real_m_traj < real_m_frame

    m = max(memory_usage(make_copy_2_frames)) - real_m_traj
    print ('make_copy_2_frames', m)

    m = max(memory_usage(perform_action_radgyr))
    print ('perform_action_radgyr', m - real_m_traj)
    assert (m - real_m_traj) < real_m_traj
    m_radgyr = m - real_m_traj

    m = max(memory_usage(perform_action_radgyr_dset_view))
    print ('perform_action_radgyr', m - real_m_traj)

    m = max(memory_usage(perform_action_radgyr_dset_copy))
    print ('perform_action_radgyr', m - real_m_traj)

