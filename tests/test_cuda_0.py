# this file was adapted from pycuda test: test_gpuarray.py)
#! /usr/bin/env python

import sys

def has_pycuda():
    try:
        import pycuda 
        return True
    except ImportError:
        return False

print (has_pycuda())

if not has_pycuda():
    print ("does not have pycuda. Quite")
    sys.exit()

import unittest # pragma no test
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

import numpy as np
import numpy.linalg as la

try:
    from pycuda.tools import mark_cuda_test
    from pycuda.characterize import has_double_support
except:
    pass


if has_pycuda():
    import pycuda.gpuarray as gpuarray
    import pycuda.driver as drv
    from pycuda.compiler import SourceModule


class TestGPUArray(unittest.TestCase):
    disabled = not has_pycuda()

    @mark_cuda_test
    def test_load_frame_coords_to_gpu(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        a = np.asarray(traj[0].coords).astype(np.float32)

        a_gpu = gpuarray.to_gpu(a)
        a_cpu = a_gpu.get()
        assert isinstance(a_gpu, gpuarray.GPUArray) == True
        assert_almost_equal(a_cpu, traj[0].coords)
        assert isinstance(a_cpu, np.ndarray) == True

if __name__ == "__main__":
    unittest.main()
