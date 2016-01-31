.. _trajectory_exercise:


Trajectory warm up
==================

.. contents::

Register to load from disk
--------------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj

Register to load from disk with frame stride
--------------------------------------------

.. ipython:: python

    # specify start=2, stop=50, step=2
    pt.iterload('tz2.nc', 'tz2.parm7', frame_slice=(2, 50, 2))

    # skip every 2 frames with stride option
    pt.iterload('tz2.nc', 'tz2.parm7', stride=2)


Register to load several files from disk
----------------------------------------

.. ipython:: python

    pt.iterload(['tz2.0.nc', 'tz2.1.nc'], 'tz2.parm7')

Register to load several files from disk with frame stride
----------------------------------------------------------

.. ipython:: python

    pt.iterload(['tz2.0.nc', 'tz2.1.nc'], 'tz2.parm7', frame_slice=[(0, 8, 2), (0, 5, 1)])

    # similiar to cpptraj style
    # parm tz2.parm7
    # trajin tz2.0.nc 1 8 2
    # trajin tz2.1.nc 1 5 1

.. ipython:: python

    # please note that pytraj skip last frame (to follow Python's convention).
    # if you specify (0, 5), pytraj will take frames from 0 to 4 (skip 5)

.. ipython:: python

    # use stride to simplify loading if do not need to specif start and stop frames
    pt.iterload(['tz2.0.nc', 'tz2.1.nc'], 'tz2.parm7', stride=3)

Write to disk
-------------

.. ipython:: python

    pt.write_traj('my_traj.nc', traj, overwrite=True)

Extract coordinates from PDB structure
--------------------------------------

.. ipython:: python

    pdb = pt.fetch_pdb('1ES7')

    # get coordinates for CA atoms
    pdb['@CA'].xyz

    # another way with more memory efficient
    pt.get_coordinates(pdb, mask='@CA')

Write PDB file to netcdf format (require Netcdf lib)
----------------------------------------------------

.. ipython:: python

    pdb = pt.fetch_pdb('1ES7')
    pt.write_traj('mypdb.nc', traj=pdb, overwrite=True)

    # can also use the shortcut
    pdb.save('mypdb.nc')


Fancy indexing
--------------

.. ipython:: python

    traj[0]
    traj[:2]
    traj[:2, '@CA']

Iterate whole trajectory
------------------------

.. ipython:: python

    for frame in traj:
        # do something with frame
        pass
    frame

Iterate a part of trajectory
----------------------------

- with stop value

.. ipython:: python

    for frame in pt.iterframe(traj, stop=5):
        print(frame)

- with given frame indices

.. ipython:: python

    for frame in pt.iterframe(traj, frame_indices=[0, 5, 20, 50]):
        print(frame)

- with given mask

.. ipython:: python

    for frame in pt.iterframe(traj, frame_indices=[0, 5, 20, 50], mask='@CA'):
        print(frame)

Perform calculation
-------------------

.. ipython:: python
    
    # rmsd to first frame with mask='@CA'
    # python starts counting from 0
    pt.rmsd(traj, ref=0, mask='@CA')

Convert data to pandas DataFrame
--------------------------------

.. ipython:: python

    df = pt.multidihedral(traj, resrange='3-7', dtype='dataframe')
    type(df)
    df.head(1)
    df.tail(1)

Convert to different file format
--------------------------------

.. ipython:: python

    # to DCD format
    pt.write_traj('traj.dcd', traj, overwrite=True)


Combine with cpptraj command style
-----------------------------------

.. ipython:: python

    pt.compute(['rms', 'radgyr @CA nomax', 'distance :3 :7'], traj)

Convert out-of-core TrajectoryIterator to in-memory Trajectory
--------------------------------------------------------------

.. ipython:: python

    traj2 = traj[:]
    # apply any transformations

    # superpose to first frame
    pt.superpose(traj2)

    # use the same syntax to perform calculation
    pt.rmsd(traj2, ref=0)

Parallel computing with TrajectoryIterator
------------------------------------------

.. ipython:: python

    # serial: pt.rmsd(traj)

    # parallel
    pt.pmap(pt.radgyr, traj, n_cores=3)

    # chain a series of cpptraj's commands
    pt.pmap(['radgyr nomax', 'molsurf @CA', 'multidihedral resrange 3-4 psi phi'], traj, n_cores=4)
