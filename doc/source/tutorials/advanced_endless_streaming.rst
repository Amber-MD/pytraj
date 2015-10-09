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
    from pytraj.datasets import DatasetList as CpptrajDatasetList
    from pytraj import ActionList

    # create pytraj.TrajectoryIterator, no data is actually loaded yet.
    traj = pt.iterload("tz2.nc", "tz2.parm7")

    # create a list of commands (cpptraj's style)
    # advantage: endless streaming calculation, even TB of data.
    commands = ['autoimage',
                'distance :3 :7',
                'distance :3 :10',
                'vector :2 :3',
                ]

    # create a CpptrajDatasetList object to hold all datasets
    dslist = CpptrajDatasetList()
    dslist

    # create an ActionList object to hold all actions
    actlist = ActionList(commands, traj.top, dslist=dslist)

    # perform the actions, only a single frame is loaded
    for frame in traj:
        actlist.do_actions(frame)
    dslist
    for d in dslist:
        print(d)
