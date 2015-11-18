Weird things in pytraj
======================

Explicit copy
-------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    xyz = [frame.xyz for frame in traj]

    # coords of 1st and last frames are the same and are [0.0, ...]
    print(xyz[0][0])
    print(xyz[-1][0])

    # need to explicitly copy
    xyz = [frame.xyz.copy() for frame in traj]

    # coords of 1st and last frames are the same
    print(xyz[0][0])
    print(xyz[-1][0])

Why does this happen? When iterating ``pytraj.TrajectoryIteatory`` (by pytraj.iterload),
pytraj create a Frame object behaving like a buffer. When getting to next snapshot, new
coordinates with box will replace the old ones. When reaching final snapshot, the Frame is
free and all coordinates are set to 0. We intend to design this behaviour since data
copying is expensive. If you want to keep the coordinates, please make sure to 'copy'.

Or you can use ``pytraj.get_coordinates``

.. ipython:: python

    xyz = pt.get_coordinates(traj, frame_indices=[0, 5, 6])
