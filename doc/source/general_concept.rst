.. _general_concept:

General concept
---------------

In most case, to use ``pytraj``, you just need to rember this format

.. code-block:: python

    import pytraj as pt
    traj = pt.iterload(...)
    data = pt.method_name(traj, mask=mask, ref=ref, frame_indices=frame_indices ...) 

* ``traj`` can be a Trajectory-like (having 3D coordinates and a Topology) or it can be
  any iterable object that produces :class:`pytraj.Frame` (a single snapshot).

* ``ref`` can be a single Frame or a Trajectory. If it's a Trajectory,
  first conformation is always picked up.

* ``mask`` follows Amber mask syntax (eg. :3-18@CA) or an atom index array (eg. [0, 3, 5]). If ``mask`` is a string (amber mask), the index is 1-based (counting from 1) and 
 if ``mass`` is an array-like, the index is 0-based (counting from 0). 

* ``frame_indices`` is frame indices for calculation. It's optional. If no ``frame_indices`` is provided, the calculation will be performed for whole trajectory

Real world example

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

.. ipython:: python
 
    # calculate radius of gyration for whole trajectory, all atoms
    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj

    data = pt.radgyr(traj)

    # calculate radius of gyration for whole trajectory, @CA atoms
    data = pt.radgyr(traj, mask='@CA')

    # calculate radius of gyration for whole trajectory, for atom indices of [0, 3, 5, 100]
    data = pt.radgyr(traj, mask=[0, 3, 5, 100])

    # calculate radius of gyration for whole trajectory, all atoms, only for 0, 7, and 8-th frames
    data = pt.radgyr(traj, frame_indices=[0, 7, 8])

    # calculate rmsd with reference as 0-th frame, all atoms
    data = pt.rmsd(traj, ref=0)

    # calculate rmsd with reference as 0-th frame, backbone heavy atoms
    data = pt.rmsd(traj, ref=0, mask='@C,CA,N,O')

    # calculate rmsd with reference as last frame, CA atoms
    data = pt.rmsd(traj, ref=-1, mask='@CA')
    data
