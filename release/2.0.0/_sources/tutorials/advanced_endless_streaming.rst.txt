.. _advanced_endless_streaming:

Endless streaming calculation
=============================

.. note:: syntax might be changed in future.

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

|

.. ipython:: python

    import pytraj as pt
    from pytraj.datasets import CpptrajDatasetList

    # create pytraj.TrajectoryIterator, no data is actually loaded yet.
    traj = pt.iterload("tz2.ortho.nc", "tz2.ortho.parm7")
    traj

    # create a list of commands (cpptraj's style)
    # advantage: endless streaming calculation, even TB of data.
    commands = ['autoimage',
                'distance :3 :7',
                'distance :3 :10',
                'vector :2 :3',
                'vector box',
                'vector ucellx',
    ]

    # create a CpptrajDatasetList object to hold all datasets
    dslist = CpptrajDatasetList()
    print(dslist)

    # create an ActionList object to hold all actions
    actlist = pt.ActionList(commands, top=traj.top, dslist=dslist)

    # perform the actions, only a single frame is loaded
    # data is saved to dslist
    for frame in traj:
        actlist.compute(frame)
        # you can plug your own function here to
        # your_funct(frame, ...)
    print(dslist)

    # get values for each Dataset
    for d in dslist: print(d)

    # get raw data for the last DatsetVector (x values of unitcells ('vector ucellx'))
    print(dslist[-1].values)
