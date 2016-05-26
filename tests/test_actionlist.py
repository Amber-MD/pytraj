#!/usr/bin/env python

from __future__ import print_function
import unittest
import pytraj as pt
import numpy as np
from pytraj import adict, allactions
from pytraj import ArgList, Trajectory, Frame
from pytraj.utils import eq, aa_eq
from pytraj.c_action import c_action as CA
from pytraj.datasets import DatasetList as CpptrajDatasetList
from pytraj.datafiles.datafiles import DataFileList
from pytraj import ActionList
from pytraj import Pipeline
from pytraj.testing import cpptraj_test_dir


class TestActionList(unittest.TestCase):

    def test_distances(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")[:]

        trajin = pt.datafiles.tc5b_trajin + """
        distance @CB @CA
        distance @CA @H
        """

        cout = pt.datafiles.load_cpptraj_output(trajin)[1:]

        mask_list = ('@CB @CA', '@CA @H')
        dslist = pt.calc_distance(traj, mask_list)
        dslist3_0 = pt.calc_distance(traj, mask_list[0])
        dslist3_1 = pt.calc_distance(traj, mask_list[1])

        # compare to cpptraj output
        aa_eq(dslist.flatten(), cout.values.flatten())
        aa_eq(dslist3_0, dslist[0])
        aa_eq(dslist3_1, dslist[1])

    def test_run_0(self):
        # load traj
        farray = pt.load(filename="./data/tz2.truncoct.nc",
                         top="./data/tz2.truncoct.parm7")[:2]
        fold = farray.copy()

        act = allactions.Action_Image()
        ptrajin = """
        center :2-11
        image center familiar com :6
        """

        # create 'strip' action
        stripact = allactions.Action_Strip()

        # creat datasetlist to hold distance data
        dsetlist = CpptrajDatasetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()

        top = farray.top

        # add two actions: Action_Strip and Action_Distance
        alist.add(allactions.Action_Center(), ArgList(":2-11"), top=top)
        alist.add(allactions.Action_Image(),
                  ArgList("center familiar com :6"),
                  top=top)

        # do checking
        alist.setup(top)

        farray2 = Trajectory()
        frame0 = Frame()
        # testing how fast to do the actions

        # loop all frames
        # use iterator to make faster loop
        # don't use "for i in range(farray.n_frames)"
        for frame in farray:
            # perform actions for each frame
            # we make a copy since we want to keep orginal Frame
            frame0 = frame.copy()
            alist.compute(frame0)

            # we need to keep the modified frame in farray2
            farray2.append(frame0)

        # make sure that Action_Strip does its job in stripping
        assert farray2.n_frames == farray.n_frames

        fsaved = pt.iterload(cpptraj_test_dir + "/Test_Image/image4.crd.save",
                             "data/tz2.truncoct.parm7")
        assert fsaved.n_frames == 2

    def test_run_1(self):
        # load traj
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        dslist = CpptrajDatasetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()
        # add two actions: Action_Dihedral and Action_Distance
        alist.add(adict['distance'],
                  ":2@CA :10@CA out ./output/_dist.out", traj.top,
                  dslist, dflist)
        alist.add(adict['dihedral'],
                  ":2@CA :3@CA :4@CA :5@CA out ./output/_dih.out",
                  traj.top, dslist, dflist)

        # using string for action 'dssp'
        alist.add('dssp', "out ./output/_dssp_alist.out", traj.top,
                  dslist, dflist)
        alist.add('matrix', "out ./output/_mat_alist.out", traj.top,
                  dslist, dflist)
        # does not work with `strip` (output traj have the same n_atoms as originl traj)
        # turn off for now
        # Error: Could not get associated topology for ./output/test_trajout.nc
        # alist.compute([traj[[0, 1]], traj, traj.iterchunk(chunksize=4,
        #                                                     stop=8),
        #                  traj.iterframe()])

    def test_run_2(self):
        # load traj
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        dslist = CpptrajDatasetList()
        dflist = DataFileList()

        # creat ActionList to hold actions
        alist = ActionList()
        alist.add(adict['distance'],
                  ":2@CA :10@CA out ./output/_dist.out", traj.top,
                  dslist, dflist)
        alist.compute([traj.iterchunk()])
        print('dslist', dslist[0].size)
        assert len(dslist) == 1
        assert dslist[0].size == traj.n_frames

    def test_run_3(self):
        dslist = CpptrajDatasetList()
        actlist = ActionList()
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mask_list = ['@CB @CA @N', '@CA @H @N']

        for mask in mask_list:
            actlist.add(CA.Action_Angle(),
                        mask,
                        traj.top,
                        dslist=dslist)
        actlist.compute(traj)

        dslist2 = pt.calc_angle(traj, mask_list)

        dslist3_0 = pt.calc_angle(traj, mask_list[0])
        dslist3_1 = pt.calc_angle(traj, mask_list[1])
        aa_eq(dslist3_0, dslist[0].to_ndarray())
        aa_eq(dslist3_1, dslist[1].to_ndarray())

        aa_eq(dslist3_0, dslist[0].to_ndarray())
        aa_eq(dslist3_1, dslist[1].to_ndarray())

    def test_run_4(self):
        dslist = CpptrajDatasetList()
        actlist = ActionList()
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mask_list = ['@CB @CA @N @H', '@CA @H @N @H=']

        for mask in mask_list:
            actlist.add(CA.Action_Dihedral(),
                        mask,
                        traj.top,
                        dslist=dslist)
        actlist.compute(traj)

        dslist2 = pt.calc_dihedral(traj, mask_list)

        dslist3_0 = pt.calc_dihedral(traj, mask_list[0])
        dslist3_1 = pt.calc_dihedral(traj, mask_list[1])
        aa_eq(dslist3_0, dslist2[0])
        aa_eq(dslist3_1, dslist2[1])

        aa_eq(dslist3_0, dslist[0].to_ndarray())
        aa_eq(dslist3_1, dslist[1].to_ndarray())

    def test_run_5(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mask_list = ('@CB @CA', '@CA @H')
        dslist = CpptrajDatasetList()
        actlist = ActionList()

        for mask in mask_list:
            actlist.add(CA.Action_Distance(),
                        mask,
                        traj.top,
                        dslist=dslist)
        actlist.compute(traj)

        dslist2 = pt.calc_distance(traj, mask_list)
        aa_eq(dslist.values, dslist2)

    def test_6(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        mask_list = ('@CB @CA', '@CA @H')
        dslist = pt.calc_distance(traj, mask_list)
        dslist3_0 = pt.calc_distance(traj, mask_list[0])
        dslist3_1 = pt.calc_distance(traj, mask_list[1])

        aa_eq(dslist3_0, dslist[0])
        aa_eq(dslist3_1, dslist[1])

    def test_constructor_from_command_list_TrajectoryIterator(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        commands = ['rmsd @CA', 'distance :3 :7', 'distance     :3 :7',
                    'vector :2 :3']

        dslist = CpptrajDatasetList()
        actlist = ActionList(commands, traj.top, dslist=dslist)

        d0 = dslist.add('ref_frame', 'my_ref')
        d0.top = traj.top
        d0.add_frame(traj[3])

        for frame in traj:
            actlist.compute(frame)

        aa_eq(pt.rmsd(traj, mask='@CA'), dslist[0])
        aa_eq(pt.distance(traj, ':3 :7'), dslist[1])
        aa_eq(pt.distance(traj, ':3 :7'), dslist[2])
        aa_eq(
            pt.vector.vector_mask(
                traj(rmsfit=(0, '@CA')),
                ':2 :3'),
            dslist[3].values)

    def test_constructor_from_command_list_Trajectory(self):
        '''mutable Trajectory'''
        # use `load` method rather `iterload`
        traj = pt.load("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        # make sure no space-sensitivity
        # make the code (commands) ugly is my intention.
        commands = [
            'autoimage ',
            'autoimage',
            'rmsd @CA',
            'distance :3 :7',
            'distance     :3 :7',
            'vector :2 :3',
            '  distance :3 :7',
            'rms @C,N,O',
        ]

        dslist = CpptrajDatasetList()
        actlist = ActionList(commands, traj.top, dslist=dslist)

        for frame in traj:
            actlist.compute(frame)

        aa_eq(pt.rmsd(traj, mask='@CA'), dslist[0])
        aa_eq(pt.distance(traj, ':3 :7'), dslist[1])
        aa_eq(pt.distance(traj, ':3 :7'), dslist[2])
        # do not need to perform rmsfit again.
        aa_eq(pt.vector.vector_mask(traj, ':2 :3'), dslist[3].values)
        aa_eq(pt.distance(traj, ':3 :7'), dslist[4])
        aa_eq(pt.rmsd(traj, mask='@C,N,O'), dslist[5])

    def test_constructor_from_command_list_TrajectoryIterator_no_DatasetList(
            self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        commands = ['rmsd @CA', 'distance :3 :7', 'distance     :3 :7',
                    'vector :2 :3']

        actlist = ActionList(commands, top=traj.top)

        for frame in traj:
            actlist.compute(frame)

        aa_eq(pt.rmsd(traj, mask='@CA'), actlist.data[0])
        aa_eq(pt.distance(traj, ':3 :7'), actlist.data[1])
        aa_eq(pt.distance(traj, ':3 :7'), actlist.data[2])
        aa_eq(
            pt.vector.vector_mask(
                traj(rmsfit=(0, '@CA')),
                ':2 :3'),
            actlist.data[3].values)

    def test_modify_frame(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        dslist = CpptrajDatasetList()
        dslist.add_new('topology', name='mytop')

        # add a new topology
        dslist[0].data = pt.strip(traj.top, ':WAT')
        commands = ['autoimage', 'strip :WAT', 'createcrd mycrd', ]

        actlist = ActionList(commands, top=traj.top, dslist=dslist)

        for frame in traj:
            actlist.compute(frame)

        aa_eq(dslist['mycrd'].xyz,
              pt.get_coordinates(traj,
                                 mask='!:WAT',
                                 autoimage=True))

    def test_modify_frame_use_Pipeline(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        dslist = CpptrajDatasetList()
        dslist.add_new('topology', name='mytop')

        # add a new topology
        dslist[0].data = pt.strip(traj.top, ':WAT')
        commands = ['autoimage', 'strip :WAT', 'createcrd mycrd', ]

        actlist = Pipeline(commands, top=traj.top, dslist=dslist)

        for frame in traj:
            actlist.compute(frame)

        aa_eq(dslist['mycrd'].xyz,
              pt.get_coordinates(traj,
                                 mask='!:WAT',
                                 autoimage=True))

    def test_combine_with_frame_iterator(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        dslist = CpptrajDatasetList()

        commands = ['autoimage', 'rms', ]

        actlist = ActionList(commands, top=traj.top, dslist=dslist)

        def get_frameiter(actlist, traj):
            for frame in traj:
                actlist.compute(frame)
                yield frame

        def do_extra(fi):
            a = []
            for frame in fi:
                frame.xyz = frame.xyz + 2.
                a.append(frame.copy())
            return a

        new_list = do_extra(get_frameiter(actlist, traj))
        t0 = traj[:].autoimage().superpose()
        t0.xyz += 2.
        aa_eq(np.array([frame.xyz for frame in new_list]), t0.xyz)

    def test_combine_cpptraj_iterating_with_pytraj(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        commands = ['autoimage', 'rms', ]

        dslist = CpptrajDatasetList()
        actlist = ActionList(commands, top=traj.top, dslist=dslist)

        def get_fi(actlist, traj):
            '''create a frame iterator with pre-processed by cpptraj
            '''
            for frame in traj:
                actlist.compute(frame)
                yield frame

        ref = traj[3]
        pt.autoimage(ref, top=traj.top)
        fi = get_fi(actlist, traj)
        rmsd_nofit_after_fitting = pt.rmsd_nofit(fi, ref=ref, top=traj.top)

        t0 = traj[:].autoimage().superpose()
        saved_rmsd_ = pt.rmsd_nofit(t0, ref=ref)
        aa_eq(rmsd_nofit_after_fitting, saved_rmsd_)

    def test_pipe(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        # from TrajectoryIterator
        fi = pt.pipe(traj, ['autoimage', 'rms'])
        xyz = np.array([frame.xyz.copy() for frame in fi])
        t0 = traj[:].autoimage().superpose()
        aa_eq(xyz, t0.xyz)

        # from FrameIterator
        fi = pt.pipe(traj(), ['autoimage', 'rms'])
        xyz = np.array([frame.xyz.copy() for frame in fi])
        t0 = traj[:].autoimage().superpose()
        aa_eq(xyz, t0.xyz)

        # from FrameIterator with indices
        fi = pt.pipe(traj(0, 8, 2), ['autoimage', 'rms'])
        xyz = np.array([frame.xyz.copy() for frame in fi])
        t0 = traj[:8:2].autoimage().superpose()
        aa_eq(xyz, t0.xyz)

        # from TrajectoryIterator, cpptraj's command style
        fi = pt.pipe(traj, '''
        autoimage
        rms''')
        xyz = np.array([frame.xyz.copy() for frame in fi])
        t0 = traj[:].autoimage().superpose()
        aa_eq(xyz, t0.xyz)

    def test_reference(self):
        traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")

        # store reference
        dslist = CpptrajDatasetList()
        ref = dslist.add_new('reference')
        ref.top = traj.top
        ref.append(traj[3])

        fi = pt.pipe(traj,
                     ['autoimage', 'rms refindex 0 @CA'],
                     dslist=dslist)
        xyz = np.array([frame.xyz.copy() for frame in fi])
        t0 = (traj[:]
              .autoimage()
              .superpose(ref=traj[3], mask='@CA'))
        print('ok')
        aa_eq(xyz, t0.xyz)

        t1 = traj[:].autoimage()
        aa_eq(pt.rmsd(t1, ref=traj[3], mask='@CA'), dslist[-1].values)


if __name__ == "__main__":
    unittest.main()
