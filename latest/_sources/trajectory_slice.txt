.. _trajectory_slice:

Fancy indexing of trajectory
============================

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj


.. contents::

use `[ ]` slicing notation to load all given frames into memory at once
-----------------------------------------------------------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    print(traj)
    
    # get 1st frame
    traj[0]

    # get last frame
    traj[-1]

    # get frame 2 to 8, skip every 2
    traj[2:8:2]

    # strip all atoms but CA
    traj['@CA']

    # take coords for 1st residues
    traj[':1']

    # get frames 1, 3, 7 with only CA atoms
    traj[[1, 3, 7], '@CA']

    # get frame 2 to 8, skip every 2 and keep only non-H atoms
    traj[2:8:2, '!@H=']

    # get frame 2 to 8, skip every 2, reverse from last to begin
    # strip water
    traj[8:2:-2, '!:WAT']

    # skip every 2 frames
    traj[::2]


use `( )` notation to create frame iterator, load only a single frame to memory if needed
-----------------------------------------------------------------------------------------

.. ipython:: python

    traj(0, 8, 2)
    for frame in traj(0, 8, 2): print(frame)

    # iterating with autoimage
    for frame in traj(0, 8, 2, autoimage=True): print(frame)

    # iterating with rmsfit to first frame
    for frame in traj(0, 8, 2, rmsfit=0): print(frame)

    # iterating all frames and keep only CA atoms
    for frame in traj(mask='@CA'): print(frame)

    # iterating all frames and strip waters
    for frame in traj(mask='!:WAT'): print(frame)
