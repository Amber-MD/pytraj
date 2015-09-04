.. _tips_iterload:

Tips for out-of-core loading
----------------------------

.. note:: when calling ``iterload``, no data is loaded to memory if not needed.

Register to load up to 50-th frame and skip every 2 frames

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7', frame_slice=(0, 50, 2))
    traj

Register to load all frames

.. ipython:: python

    import pytraj as ptraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj

    # since we register to load all frames, if we just need up to 50-th frame
    # we juse ( ) notation
    for frame in traj(0, 50): pass
    print(frame)
