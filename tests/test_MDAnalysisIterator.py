from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca
from pytraj.compat import zip

class Test(unittest.TestCase):
    @test_if_having("MDAnalysis")
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from MDAnalysis import Universe
        from pytraj.trajs.TrajectoryMDAnalysisIterator import (
                TrajectoryMDAnalysisIterator as MDIterator)

        u = Universe(traj.top.filename, traj.filename, format='mdcrd', topology_format='prmtop')
        traj_converted = mdio.load_MDAnalysis(u, top=traj.top)
        aa_eq(traj.xyz, traj_converted.xyz)
        aa_eq(traj[0].xyz, traj_converted[0].xyz)
        print (u.atoms)
        u_traj = u.trajectory

        titer = MDIterator(u, top=traj.top)
        titer2 = mdio.load_MDAnalysisIterator(u)
        aa_eq(traj.xyz, titer.xyz)
        aa_eq(traj.xyz, titer2.xyz)
        print (titer)
        assert titer.n_atoms == u_traj.numatoms
        assert titer.n_frames == u_traj.numframes

        # make sure titer.top is Topology object
        assert isinstance(titer.top, Topology)
        print (titer[0][0])

        # make sure we get correct frame with given index
        aa_eq(traj_converted[0].xyz, titer[0].xyz)
        aa_eq(traj_converted[0].xyz, titer2[0].xyz)

        # test slicing
        aa_eq(titer[0:10:2].xyz, traj[0:10:2].xyz)
        aa_eq(titer2[0:10:2].xyz, traj[0:10:2].xyz)

        # make sure we can reprodue pytraj' xyz coords
        for f0, f1 in zip(titer, traj):
            #print (f0[0])
            assert isinstance(f0, Frame) == True
            aa_eq(f0.xyz, f1.xyz)

        # let's do some analysis
        d_mda = pyca.search_hbonds(titer, dtype='ndarray')
        d_mda2 = pyca.search_hbonds(titer2, dtype='ndarray')
        d_traj = pyca.search_hbonds(traj, dtype='ndarray')
        aa_eq(d_mda, d_traj)

        # FIXME : failed assertion
        #aa_eq(d_mda2, d_traj)
        print (d_mda2)
        print (d_traj)

    @no_test
    @test_if_having("MDAnalysis")
    def test_1(self):
        # FIXME: I don't know how MDAnalysis rewind this DCD file
        # DCD and PSF: always got trouble with MDAnalysis
        # for those files
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import DCD, PSF
        u = Universe(PSF, DCD)
        t = mdio.load_MDAnalysisIterator(u)
        traj = mdio.load(DCD, PSF)
        t_in_mem = mdio.load_MDAnalysis(u)
        aa_eq(t.xyz, traj.xyz)

        # assert failed
        # IOError: Error reading frame from DCD file
        #aa_eq(t[:2].xyz, traj[:2].xyz)
        d0 = traj.search_hbonds()
        print (d0)
        # make sure no segfault
        d1 = t.search_hbonds()
        d1 = t.search_hbonds()
        d1 = t.search_hbonds()
        d1 = t.search_hbonds()
        d1 = t.search_hbonds()
        print (d1.keys())
        # FIXME: how MDA handle DCD file?
        # assertion failed.
        #aa_eq(d0.to_ndarray(), d1.to_ndarray())

        ## try another action: COM
        d0 = traj.calc_COM().to_ndarray()
        d1 = t.calc_COM().to_ndarray()
        aa_eq(d0, d1)

        # try another action: COG
        d0 = traj.calc_COG().to_ndarray()
        d1 = t.calc_COG().to_ndarray()
        aa_eq(d0, d1)

        # try another action: RMSD
        d0 = pyca.calc_rmsd(traj, ref=traj[-1])
        # FIXME: segfault
        #d1 = pyca.calc_rmsd(t, ref=traj[-1])
        #aa_eq(d0, d1)

        # try another action: DSSP: segfault
        # FIXME: segfault
        #d0 = pyca.calc_dssp(traj, dtype='ndarray')
        #d1 = pyca.calc_dssp(t, dtype='ndarray')
        #d2 = pyca.calc_dssp(t_in_mem, dtype='ndarray')
        #aa_eq(d0, d1)
        #aa_eq(d0, d2)

        # try another action: rms2d
        #d0 = pyca.calc_pairwise_rmsd(traj, dtype='ndarray')
        # FIXME: segfault
        #d1 = pyca.calc_pairwise_rmsd(t, dtype='ndarray')
        # FIXME: segfault
        #d2 = pyca.calc_pairwise_rmsd(t_in_mem, dtype='ndarray')
        #aa_eq(d0, d1)
        #aa_eq(d0, d2)

        # just try not to mak segfault
        # FIXME: segfault
        #pyca.do_clustering(t, "kmeans clusters 5 @CA")
        #pyca.do_clustering(t[:], "kmeans clusters 5 @CA")
        #pyca.do_clustering(t_in_mem, "kmeans clusters 5 @CA")
        #print (t)

    @test_if_having("MDAnalysis")
    def test_2(self):
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import XYZ_psf, XYZ_bz2
        u = Universe(XYZ_psf, XYZ_bz2)
        print (u)
        t = mdio.load_MDAnalysisIterator(u)
        traj = mdio.load_MDAnalysis(u)
        aa_eq(t.xyz, traj.xyz)
        aa_eq(t[:2].xyz, traj[:2].xyz)
        fa0 = t[:]
        fa1 = traj[:]
        aa_eq(fa0['@CA'].xyz, fa1['@CA'].xyz)

        # try another action: rms2d
        d0 = pyca.calc_pairwise_rmsd(traj, dtype='ndarray')
        d1 = pyca.calc_pairwise_rmsd(t, dtype='ndarray')
        aa_eq(d0, d1)

        # just try not to mak segfault
        pyca.do_clustering(t, "kmeans clusters 5 @CA")
        pyca.do_clustering(t[:], "kmeans clusters 5 @CA")

        # try another action: RMSD
        d0 = pyca.calc_rmsd(traj, ref=traj[-1])
        d1 = pyca.calc_rmsd(t, ref=traj[-1])
        aa_eq(d0, d1)

    @test_if_having("MDAnalysis")
    def test_3(self):
        # GRO
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import GRO, TRR
        u = Universe(GRO, TRR)

        print (u)
        t = mdio.load_MDAnalysisIterator(u)
        print (t)
        traj = mdio.load_MDAnalysis(u)
        print (traj)
        xyz = t.xyz
        aa_eq(t.xyz, traj.xyz)
        fa0 = t[:]
        fa1 = traj[:]
        aa_eq(fa0['@CA'].xyz, fa1['@CA'].xyz)

        # check if segfault
        t[:2]
        t[:2]
        t[:2]
        t[:2]
        aa_eq(t[:2].xyz, traj[:2].xyz)
        pyca.calc_dssp(t)
        pyca.calc_rmsd(t, ref=traj[0])


if __name__ == "__main__":
    unittest.main()
