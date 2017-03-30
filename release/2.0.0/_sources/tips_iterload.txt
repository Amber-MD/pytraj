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

    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj
    # since we registered to load all frames, if we only need up to 50-th frame for our
    # calculation, we can use ( ) notation
    for frame in traj(0, 50): pass
    # print the last frame
    print(frame)
