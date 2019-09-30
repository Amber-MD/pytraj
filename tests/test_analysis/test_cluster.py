from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from pytraj import cluster
from pytraj.utils.context import capture_stdout

# local
from utils import fn

tz2_trajin = fn('tz2.nc')
tz2_top = fn('tz2.parm7')


def test_ClusteringDataset():
    traj = pt.load(tz2_trajin, tz2_top)
    x = pt.cluster.kmeans(traj, n_clusters=5, metric='rms', mask='@CA')
    assert x.n_frames == 101
    assert list(x.centroids) == [24, 101, 76, 13, 9]
    aa_eq(x.fraction, [0.485, 0.238, 0.139, 0.079, 0.059], decimal=3)


def test_cluster_kmeans(tmpdir):

    kmeans = cluster.kmeans
    traj = pt.iterload(tz2_trajin, tz2_top)

    cm = 'cluster kmeans clusters 10 @CA rms randompoint kseed 2 {0} noinfo'

    for sieve in [2, 3, 4]:
        sieve_str = "sieve {0} sieveseed 12".format(sieve)
        command = cm.format(sieve_str) + ' cpopvtime cpptraj_cpopvtime.agr normframe'
        state = pt.load_cpptraj_state(command, traj)
        state.run()

        data = kmeans(
            traj,
            n_clusters=10,
            kseed=2,
            random_point=True,
            metric='rms',
            mask='@CA',
            options=sieve_str + ' cpopvtime pytraj_cpopvtime.agr normframe')
        assert data.n_frames == traj.n_frames
        assert os.path.exists('cpptraj_cpopvtime.agr')
        # Make sure pytraj could write data file.
        assert os.path.exists('pytraj_cpopvtime.agr')


def test_cluster_dbscan():
    command = """
    parm {}
    trajin {}
    createcrd crd1
    cluster crdset crd1 C0 @CA dbscan epsilon 1.7 minpoints 5
    """.format(tz2_top, tz2_trajin)

    with tempfolder():
        state = pt.load_cpptraj_state(command)
        with capture_stdout() as (out, _):
            state.run()
        traj = pt.iterload(tz2_trajin, tz2_top)
        data = pt.cluster.dbscan(
            traj, mask='@CA', options='epsilon 1.7 minpoints 5')
        aa_eq(state.data[-2], data.cluster_index)


def test_cluster_hieragglo():
    command = """
    parm {}
    trajin {}
    createcrd crd1
    cluster crdset crd1 C0 !@H hieragglo epsilon 0.8 averagelinkage
    """.format(tz2_top, tz2_trajin)

    with tempfolder():
        state = pt.load_cpptraj_state(command)
        with capture_stdout() as (cpp_out, _):
            state.run()
        traj = pt.iterload(tz2_trajin, tz2_top)
        data = pt.cluster.hieragglo(
            traj, mask='!@H=', options='epsilon 0.8 averagelinkage')
        aa_eq(state.data[-2], data.cluster_index)
