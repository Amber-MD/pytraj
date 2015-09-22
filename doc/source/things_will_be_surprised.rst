Weird things in pytraj
======================

Explit copy
-----------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    xyz = [frame.xyz for frame in traj]

    # coords of 1st and last frames are the same
    print(xyz[0][0])
    print(xyz[-1][0])

    # need to explicitly copy
    xyz = [frame.xyz.copy() for frame in traj]

    # coords of 1st and last frames are the same
    print(xyz[0][0])
    print(xyz[-1][0])
