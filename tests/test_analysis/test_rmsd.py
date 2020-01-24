#!/usr/bin/env python

import unittest
import numpy as np
import pytraj as pt

from pytraj.utils import has_
from pytraj.testing import aa_eq, cpptraj_test_dir, tempfolder
from pytraj import Trajectory
from pytraj.datasets import CpptrajDatasetList
from pytraj import AtomMask

from utils import fn, tc5b_trajin, tc5b_top, tz2_trajin, tz2_top


class TestSimpleRMSD(unittest.TestCase):
    def setUp(self):
        self.traj = pt.iterload(tc5b_trajin, tc5b_top)

    def test_fit_and_then_nofit(self):
        traj = pt.iterload(tc5b_trajin, tc5b_top)
        t0 = traj[:]
        pt.superpose(t0, ref=traj[3], mask='@CA')
        rmsd_0 = pt.rmsd_nofit(traj, ref=traj[3], mask='@CB')
        rmsd_1 = pt.rmsd(traj, ref=traj[3], mask='@CB', nofit=True)
        aa_eq(rmsd_1, rmsd_0)

    def test_rmsd_with_mask(self):
        TRAJ = pt.iterload(filename=tc5b_trajin, top=tc5b_top)
        cpptraj_rmsd = np.loadtxt(
            fn("rmsd_to_firstFrame_CA_allres.Tc5b.dat"),
            skiprows=1).transpose()[1]
        f0 = TRAJ[0]
        arr0 = np.zeros(TRAJ.n_frames)
        arr1 = np.zeros(TRAJ.n_frames)
        mask = "@CA"
        atm = AtomMask(mask)
        TRAJ.top._set_integer_mask(atm)

        for i, frame in enumerate(TRAJ):
            arr0[i] = frame.rmsd(f0, mask=mask, top=TRAJ.top)
            arr1[i] = frame.rmsd(f0, atommask=atm)

        arr2 = pt.rmsd(TRAJ, mask=mask, ref=f0)
        arr3 = pt.rmsd(TRAJ, mask=mask, ref=0)
        aa_eq(arr0, cpptraj_rmsd, decimal=3)
        aa_eq(arr1, cpptraj_rmsd, decimal=3)
        aa_eq(arr2, cpptraj_rmsd, decimal=3)
        aa_eq(arr3, cpptraj_rmsd, decimal=3)

    def testsuperpose_alias(self):
        '''testsuperpose_alias'''
        t0 = self.traj[:]
        t1 = self.traj[:]
        pt.transform(t0, ['superpose'])
        pt.transform(t1, ['rms'])
        aa_eq(t0.xyz, t1.xyz)

    def test_reference_with_different_topology_basic(self):
        traj1 = pt.iterload(filename=tc5b_trajin, top=tc5b_top)
        traj2 = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        # re-establish ActionList
        dslist = CpptrajDatasetList()
        dslist.add('reference', name='myref')

        dslist[0].top = traj2.top
        dslist[0].add_frame(traj2[0])

        actlist = pt.ActionList(
            ['rmsd @1-11 @CB ref myref'], top=traj1.top, dslist=dslist)
        for frame in traj1:
            actlist.compute(frame)

        # raise if ref_mask is given but not mask
        self.assertRaises(ValueError,
                          lambda: pt.rmsd(traj1, ref=3, ref_mask='@CB'))
        self.assertRaises(
            ValueError, lambda: pt.rmsd(traj1, ref=traj2[:1], ref_mask='@CB'))

        # assert to cpptraj
        tc5b_traj = traj1[:]
        tz2_traj = traj2[:1]

        cm = '''
        parm  {} [tc5b]
        trajin {}
        parm {} [tz2]
        reference {} parm [tz2] 1 [myref]
        rms myrmsd ref [myref] @1-10 @11-20
        '''.format(tc5b_top, tc5b_trajin, tz2_top, tz2_trajin)
        print(cm)
        state = pt.load_cpptraj_state(cm)
        with tempfolder():
            state.run()

        expected_rmsd = state.data[-1].values
        rmsd_data = pt.rmsd(
            tc5b_traj, mask='@1-10', ref=tz2_traj, ref_mask='@11-20')
        aa_eq(expected_rmsd, rmsd_data)

    @unittest.skipIf(not has_('mdtraj'), 'does not have mdtraj')
    def test_ComparetoMDtraj(self):
        import mdtraj as md
        traj = pt.load(filename=tc5b_trajin, top=tc5b_top)
        m_top = md.load_prmtop(tc5b_top)
        m_traj = md.load_mdcrd(tc5b_trajin, m_top)
        m_traj.xyz = m_traj.xyz * 10  # convert `nm` to `Angstrom` unit

        arr0 = pt.rmsd(traj, ref=0)
        arr1 = pt.rmsd(traj, ref=0)
        arr2 = pt.rmsd(traj, )
        a_md0 = md.rmsd(m_traj, m_traj, 0)
        aa_eq(arr0, arr1)
        aa_eq(arr0, arr2)
        aa_eq(arr0, a_md0)

        arr0 = pt.rmsd(traj, ref=-1)
        arr1 = pt.rmsd(traj, ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1)
        aa_eq(arr0, arr1)
        aa_eq(arr0, a_md)

        mask = ":3-18@CA,C"
        atm = traj.top(mask)
        arr0 = pt.rmsd(traj, ref=-1, mask=mask)
        arr1 = pt.rmsd(traj, mask=atm.indices, ref=-1)
        arr2 = pt.rmsd(traj, mask=list(atm.indices), ref=-1)
        arr3 = pt.rmsd(traj, mask=tuple(atm.indices), ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1, atm.indices)
        aa_eq(arr0, a_md)
        aa_eq(arr1, a_md)
        aa_eq(arr2, a_md)
        aa_eq(arr3, a_md)

        fa = Trajectory(traj)
        arr0 = pt.rmsd(fa, ref=-1, mask=mask)
        arr1 = pt.rmsd(fa, mask=atm.indices, ref=-1)
        arr2 = pt.rmsd(fa, mask=list(atm.indices), ref=-1)
        arr3 = pt.rmsd(fa, mask=tuple(atm.indices), ref=-1)
        a_md = md.rmsd(m_traj, m_traj, -1, atm.indices)
        aa_eq(arr0, a_md)
        aa_eq(arr1, a_md)
        aa_eq(arr2, a_md)
        aa_eq(arr3, a_md)

        fa = Trajectory(traj)
        mask = "!@H="
        atm = fa.top(mask)
        arr0 = pt.rmsd(fa, ref=4, mask=mask)
        a_md = md.rmsd(m_traj, m_traj, 4, atm.indices)

        # exclude 0-th frame for ref
        aa_eq(arr0, a_md)

    def test_list_of_masks(self):
        traj = self.traj.copy()
        mask = ['@CA', '@CB', ':3-18@CA,C']
        arr = pt.rmsd(traj, mask=mask)
        for idx, m in enumerate(mask):
            aa_eq(arr[idx], pt.rmsd(traj, mask=m))
            aa_eq(arr[idx], pt.rmsd(traj, mask=traj.top.select(m)))

        mask = ['@CA', '@CB', ':3-18@CA,C', [0, 3, 5]]
        self.assertRaises(TypeError, lambda: pt.rmsd(traj, mask=mask))

        mask_2 = [[0, 3, 6], range(50)]
        aa_eq(pt.rmsd(traj, mask=mask_2)[0], pt.rmsd(traj, mask=mask_2[0]))
        aa_eq(pt.rmsd(traj, mask=mask_2)[1], pt.rmsd(traj, mask=mask_2[1]))

        ca = pt.select('@CA', traj.top)
        cb = pt.select('@CB', traj.top)
        aa_eq(pt.rmsd(traj, mask=ca), pt.rmsd(traj, mask=[ca, cb])[0])
        aa_eq(pt.rmsd(traj, mask=cb), pt.rmsd(traj, mask=[ca, cb])[1])

    def test_raise_savematrices_if_not_dataset(self):
        traj = self.traj.copy()
        self.assertRaises(
            ValueError,
            lambda: pt.rmsd(traj, mask='@CA savematrices', dtype='ndarray'))

    def test_not_update_coordinates(self):
        traj = self.traj[:]
        data = pt.rmsd(traj, ref=3, update_coordinate=False)

        # make sure coordinates are not updated
        aa_eq(traj.xyz, self.traj.xyz)

        # make sure give the same rmsd values
        aa_eq(pt.rmsd(traj, ref=3), data)

    def test_combine_nofit_mass_nomod(self):
        cm = '''
        parm {}
        trajin {}
        rms @CA nofit mass nomod
        '''.format(fn('tz2.parm7'), fn('tz2.nc'))
        state = pt.load_cpptraj_state(cm)
        state.run()

        unmut_traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))
        mut_traj = unmut_traj[:]

        data = pt.rmsd(
            mut_traj,
            mask='@CA',
            mass=True,
            nofit=True,
            update_coordinate=False)
        aa_eq(data, state.data[-1])


class TestRMSDPerRes(unittest.TestCase):
    def test_noreference(self):
        from pytraj.datafiles import load_cpptraj_output, tz2_ortho_trajin
        traj = pt.iterload(fn("tz2.ortho.nc"), fn("tz2.ortho.parm7"))
        cout = load_cpptraj_output("""
        parm {}
        trajin {}
        rmsd first @CA perres range 2-7""".format(
            fn('tz2.ortho.parm7'), fn('tz2.ortho.nc')))
        d = pt.rmsd_perres(
            traj, ref=0, mask='@CA', resrange='2-7', dtype='ndarray')
        aa_eq(cout[1:].values, d)

    def test_reference(self):
        traj = pt.iterload(fn("tz2.truncoct.nc"), fn("tz2.truncoct.parm7"))
        txt = '''
        reference {} 2 2
        rmsd :2-11 refindex 0 perres perresout center.agr range 1 perrescenter
        '''.format(fn('tz2.truncoct.nc'))

        with tempfolder():
            state = pt.load_batch(traj, txt).run()
        # state.data has 3 datasets: ref, rmsd, rmsd perres

        # cpptraj use 2nd reference
        rmsd0 = pt.rmsd(traj, ref=1, mask=':2-11')
        rmsdperres = pt.rmsd_perres(
            traj,
            ref=1,
            mask=':2-11',
            perres_mask='*',
            resrange='1',
            perres_center=True)
        aa_eq(rmsd0, state.data[2])
        aa_eq(rmsdperres[1], state.data[3].values)

    def test_frame_indices(self):
        traj = pt.iterload(fn("tz2.truncoct.nc"), fn("tz2.truncoct.parm7"))
        traj2 = pt.iterload(
            fn("tz2.truncoct.nc"),
            fn("tz2.truncoct.parm7"),
            frame_slice=(2, 8))

        txt = '''
        reference {} 2 2
        rmsd :2-11 refindex 0 perres perresout center.agr range 1 perrescenter
        '''.format(fn('tz2.truncoct.nc'))

        state = pt.load_batch(traj2, txt)
        with tempfolder():
            state.run()

        frame_indices = range(2, 8)
        rmsd0 = pt.rmsd(traj, ref=1, mask=':2-11', frame_indices=frame_indices)
        rmsdperres = pt.rmsd_perres(
            traj,
            ref=1,
            mask=':2-11',
            perres_mask='*',
            resrange='1',
            perres_center=True,
            frame_indices=frame_indices)
        aa_eq(rmsd0, state.data[2])
        aa_eq(rmsdperres[1], state.data[3].values)


class TestRMSDnofit(unittest.TestCase):
    def test_0(self):
        tz2_ortho_trajin = fn('tz2.ortho.nc')
        tz2_ortho_top = fn('tz2.ortho.parm7')
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)

        cout = pt.datafiles.load_cpptraj_output("""
        parm {}
        trajin {}
        rms first nofit
        rms first mass
        """.format(tz2_ortho_top, tz2_ortho_trajin))
        aa_eq(pt.rmsd(traj, nofit=True), cout[1])
        aa_eq(pt.rmsd(traj, mass=True), cout[2])


class TestPairwiseRMSD(unittest.TestCase):
    def testTwoTrajTypes(self):
        '''test different metrics with different traj objects
        '''
        funclist = [pt.iterload, pt.load]
        txt = '''
        parm {}
        trajin {}
        rms2d @CA metric_holder rmsout tmp.out
        '''.format(tc5b_top, tc5b_trajin)

        for func in funclist:
            traj = func(tc5b_trajin, tc5b_top)
            for metric in ['rms', 'nofit', 'dme']:
                d0 = pt.pairwise_rmsd(traj(mask='@CA'), metric=metric)
                d1 = pt.pairwise_rmsd(traj, mask='@CA', metric=metric)
                d2 = pt.pairwise_rmsd(traj(), mask='@CA', metric=metric)

                txt0 = txt.replace('metric_holder', metric)
                state = pt.load_cpptraj_state(txt0)
                with tempfolder():
                    state.run()
                d3 = state.data[-1].values

                aa_eq(d0, d1)
                aa_eq(d0, d2)
                aa_eq(d0, d3)


class TestActionListRMSD(unittest.TestCase):
    def test_actionlist(self):
        traj = pt.iterload(tc5b_trajin, tc5b_top)
        standard_rmsd = pt.rmsd(traj, mask='@CA')

        def test_rmsd(input_traj):
            from pytraj.analysis.c_action.c_action import Action_Rmsd
            from pytraj.datasets import DatasetList
            dslist = DatasetList()
            act = Action_Rmsd()
            act.read_input('first @CA', top=input_traj.top, dslist=dslist)
            act.setup(input_traj.top)

            for frame in input_traj:
                act.compute(frame)
            return (dslist.values)

        def test_rmsd_actlist(input_traj):
            from pytraj.analysis.c_action.c_action import Action_Rmsd
            from pytraj import ActionList
            from pytraj.datasets import DatasetList

            alist = ActionList()
            dslist = DatasetList()
            act = Action_Rmsd()
            alist.add(act, 'first @CA', top=input_traj.top, dslist=dslist)

            for frame in input_traj:
                alist.compute(frame)
            return (dslist.values)

        rmsd0 = test_rmsd(traj)
        rmsd1 = test_rmsd(traj[:])
        rmsd2 = test_rmsd_actlist(traj)
        rmsd3 = test_rmsd_actlist(traj[:])
        t0 = traj[:]
        rmsd4 = test_rmsd_actlist(t0)
        aa_eq(standard_rmsd, rmsd0)
        aa_eq(standard_rmsd, rmsd1)
        aa_eq(standard_rmsd, rmsd2)
        aa_eq(standard_rmsd, rmsd3)
        aa_eq(standard_rmsd, rmsd4)


class TestSymmRmsd(unittest.TestCase):
    def test_symmrmsd(self):
        tyr_trajin = cpptraj_test_dir + '/Test_SymmRmsd/TYR.nc'
        tn = cpptraj_test_dir + '/Test_SymmRmsd/TYR.parm7'
        saved_traj = pt.iterload(
            cpptraj_test_dir + '/Test_SymmRmsd/TYR.remap.crd.save', tn)

        traj_on_disk = pt.iterload(tyr_trajin, tn)
        traj_on_mem = pt.load(tyr_trajin, tn)

        aa_eq(traj_on_disk.xyz, traj_on_mem.xyz)

        data = pt.symmrmsd(traj_on_mem, remap=True)

        cm = """
        parm {}
        trajin {}
        symmrmsd first remap myrmsd
        createcrd mycrd
        """.format(tn, tyr_trajin)

        state = pt.load_cpptraj_state(cm)
        with tempfolder():
            state.run()

        # rmsd
        aa_eq(state.data['myrmsd'].values, data)

        # coordinates
        aa_eq(state.data['mycrd'].xyz, saved_traj.xyz, decimal=3)
        aa_eq(state.data['mycrd'].xyz, traj_on_mem.xyz, decimal=3)


def test_distance_rmsd():
    traj = pt.iterload(tz2_trajin, tz2_top)
    txt = '''
    parm {}
    trajin {}
    drmsd drms_nofit out drmsd.dat
    rms rms_nofit out drmsd.dat nofit
    rms rms_fit out drmsd.dat
    drmsd drms_fit out drmsd.dat
    '''.format(tz2_top, tz2_trajin)

    state = pt.load_cpptraj_state(txt)
    with tempfolder():
        state.run()
    cpp_data = state.data[1:]

    # distance_rmsd
    data_drmsd = pt.distance_rmsd(traj)
    aa_eq(data_drmsd, cpp_data[0])
    aa_eq(pt.drmsd(traj), cpp_data[0])

    # rms_nofit
    aa_eq(cpp_data[1], pt.rmsd(traj, nofit=True))

    # rms_fit
    aa_eq(cpp_data[2], pt.rmsd(traj, nofit=False))

    # drmsd with rmsfit
    aa_eq(cpp_data[3], pt.distance_rmsd(traj(rmsfit=0), ref=traj[0]))
