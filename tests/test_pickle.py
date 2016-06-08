#!/usr/bin/env python
from __future__ import print_function
import unittest
import numpy as np
from functools import partial
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import assert_equal_topology
from pytraj.compat import zip
from pytraj.core import Box


class TestBuildAndPickleTopology(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

    def test_convert_to_dict_and_rebuild(self):
        '''test_convert_to_dict_and_rebuild
        '''
        top = self.traj.top
        d = top.to_dict()

        new_top = pt.Topology()

        MOLNUM = 0

        for idx, (aname, atype, charge, mass, resid, resname,
                  mol_number) in enumerate(zip(d['atom_name'], d[
                      'atom_type'], d['atom_charge'], d['atom_mass'], d[
                          'resid'], d['resname'], d['mol_number'])):
            atom = pt.core.Atom(name=aname,
                                type=atype,
                                charge=charge,
                                mass=mass,
                                resid=resid)
            atom.set_mol(mol_number)
            residue = pt.core.Residue(resname, resid)
            if idx == 0:
                new_top.start_new_mol()
            if mol_number > MOLNUM:
                new_top.start_new_mol()
                MOLNUM += 1
            new_top.add_atom(atom, residue)

        new_top.add_bonds(d['bond_index'])
        new_top.add_dihedrals(d['dihedral_index'])
        new_top.box = Box(top.box.values)

        assert_equal_topology(top, new_top, self.traj)

    def test_picle(self):
        '''test_picle
        '''
        pt.to_pickle(self.traj.top, 'output/new_top.pk')
        new_top = pt.read_pickle('output/new_top.pk')
        assert_equal_topology(self.traj.top, new_top, self.traj)

    def test_to_and_from_dict(self):
        cls = self.traj.top.__class__
        top = self.traj.top
        assert_equal_topology(top, cls.from_dict(top.to_dict()), self.traj)


class TestPickleFrame(unittest.TestCase):

    def test_set_mass_correctly(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        f0 = traj[0]

        f1 = f0.__class__(f0.n_atoms)
        f1.xyz[:] = f0.xyz
        assert pt.tools.rmsd(f1.mass,
                             f0.mass) > 1.0, 'must have different mass'

        f1._set_mass_from_array(f0.mass)
        aa_eq(f1.mass, f0.mass)

    def test_pickle_frame(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        f0 = traj[0]
        print(f0.mass)

        fn = 'output/f.pk'
        pt.to_pickle(f0, fn)
        f1 = pt.read_pickle(fn)
        print(f1.mass)


class TestPickleTrajectoryIterator(unittest.TestCase):

    def test_trajiter_normal(self):
        for _pickle_topology in [True, False]:
            for frame_slice in [(0, 8, 2), (0, 10, 1)]:
                traj = pt.iterload("data/Tc5b.x",
                                   "data/Tc5b.top",
                                   frame_slice=frame_slice)
                traj._pickle_topology = _pickle_topology
                pt.io.to_pickle(traj, 'output/test0.pk')
                t0 = pt.io.read_pickle('output/test0.pk')

                aa_eq(traj.xyz, t0.xyz)
                assert_equal_topology(traj.top, t0.top, traj)

    def test_trajiter_with_actionlist(self):
       traj = pt.iterload("data/tz2.ortho.nc",
                          "data/tz2.ortho.parm7")
       traj.autoimage().center('origin').superpose('@CA')
       fn = 'output/test.pk'
       pt.to_pickle(traj, fn)
       traj2 = pt.read_pickle(fn)
       print(traj2._transform_commands)
       aa_eq(traj.xyz, traj2.xyz)


def worker(rank, frame, traj):
    pt.nastruct(traj, ref=frame)


class TestPickleFrame(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("./data/Test_NAstruct/x3dna/rna.pdb")

    def test_frame(self):
        traj = pt.iterload("data/Tc5b.x", "data/Tc5b.top")
        frame = traj[0]

        pt.to_pickle(frame, 'output/frame.pk')
        frame0 = pt.read_pickle('output/frame.pk')
        aa_eq(frame0.xyz, frame.xyz)

        # test list of frames
        fname = 'output/flist.pk'
        pt.to_pickle([traj[0], traj[1]], fname)
        f01 = pt.read_pickle(fname)
        aa_eq(np.array([f.xyz for f in f01]), traj[[0, 1]].xyz)

    def test_multiprocessing(self):
        frame0 = self.traj[0]

        from multiprocessing import Pool
        func = partial(worker, frame=frame0, traj=self.traj)
        Pool(2).map(func, range(2))


class TestPickleDatasetList(unittest.TestCase):

    def test_pickle_datasetlist(self):
        traj = pt.iterload("data/Tc5b.x", "data/Tc5b.top")
        dslist = pt.multidihedral(traj)
        pt.to_pickle(dslist, 'output/ds.pk')
        dslist2 = pt.read_pickle('output/ds.pk')
        aa_eq(dslist.values, dslist2.values)


if __name__ == "__main__":
    unittest.main()
