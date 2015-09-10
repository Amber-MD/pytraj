.. _tips:

Tips
====

.. contents::

.. _process_many_files:

process many files
------------------

Normally, user needs to follow below code to process many files

.. ipython:: python

    import pytraj as pt
    import numpy as np
    template = 'tz2.%s.nc'
    flist = [template % str(i) for i in range(4)]
    print(flist)
    trajlist = [pt.iterload(fname, 'tz2.parm7') for fname in flist]

    # calculate radgyr
    data = []
    for traj in trajlist:
        data.append(pt.radgyr(traj))
    data = np.array(data).flatten()
    print(data)

However, ``pytraj`` offer shorter (and easy?) way to do

.. ipython:: python
    
    # load all files into a single TrajectoryIterator
    traj = pt.iterload('./tz2.*.nc', 'tz2.parm7')
    # perform calculation
    data = pt.radgyr(traj)
    # that's it
    print(data)

``pytraj`` will auto-allocate ``data`` too

.. ipython:: python
    
    print([t.filename.split('/')[-1] for t in trajlist])
    data = pt.radgyr(trajlist, top=trajlist[0].top)
    print(data)
    

memory saving
-------------

If memory is critical, do not load all frames into memory.

.. ipython:: python

    # do this (only a single frame will be loaded to memory)
    pt.radgyr(traj(frame_indices=[0, 200, 300, 301]))

    # rather doing (all 4 frames will be loaded to memory)
    pt.radgyr(traj[[0, 200, 300, 301]])

    traj(frame_indices=[0, 200, 300, 301])
    traj[[0, 200, 300, 301]]

See also: :ref:`trajectory_slice`

convert trajectory
------------------

.. ipython::
    
    pt.iterload('traj.nc', 'prmtop').save('traj.dcd')
