from __future__ import print_function
#import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
from pytraj.testing import make_fake_traj
import pytraj.common_actions as pyca
from memory_profiler import memory_usage
from pytraj.compat import range
from pytraj import Trajectory
import mdtraj as md
import time

fname = "./data/nogit/tip3p/md.trj"
topname = "./data/nogit/tip3p/tc5bwat.top"


def pytraj_load():
    io.load(fname, topname)


def pytraj_constructor():
    traj = Trajectory(fname, topname)


def pytraj_api():
    from pytraj import api
    traj = api.Trajectory(fname, topname)


def mdtraj_load():
    top = md.load_prmtop(topname)
    md.load_mdcrd(fname, top=top)

if __name__ == '__main__':
    import numpy as np
    from numpy import max

    m_pytraj_c = max(memory_usage(pytraj_constructor))
    print("pytraj_constructor", m_pytraj_c)
    m_pytraj_l = max(memory_usage(pytraj_load))
    print("pytraj_load", m_pytraj_l)
    m_pytraj_api = max(memory_usage(pytraj_api))
    print("pytraj_api", m_pytraj_api)
    m_mdtraj = max(memory_usage(mdtraj_load))
    print("mdtraj_load", m_mdtraj)
