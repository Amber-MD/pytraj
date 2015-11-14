.. _trajectory_excercise:


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

    # start=0, stop=8, step=2
    pt.iterload('tz2.nc', 'tz2.parm7', frame_slice=(0, 8, 2))


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

    # please note that pytraj skip last frame (to follow Python's convention).
    # if you specify (0, 5), pytraj will take frames from 0 to 4 (skip 5)


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

Do computing
------------

.. ipython:: python
    
    # rmsd to first frame with mask='@CA'
    # python starts counting from 0
    pt.rmsd(traj, ref=0, mask='@CA')

Convert data to pandas DataFrame
--------------------------------

.. ipython:: python

    df = pt.multidihedral(traj, resrange='3-7', dtype='dataframe')
    type(df)
    df.head()
    df.tail()

Convert to different file format
--------------------------------

.. ipython:: python

    # to DCD format
    pt.write_traj('traj.dcd', traj, overwrite=True)


Combine with cpptraj commmand style
-----------------------------------

.. ipython:: python

    pt.do(['rms', 'radgyr @CA', 'distance :3 :7'], traj)
