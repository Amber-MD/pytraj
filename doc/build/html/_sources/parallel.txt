.. _parallel:

Parallel calculation
====================

``pytraj`` support different types of parallel computing.

.. note:: ``pytraj`` use `python multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_, so users don't need to install extra package.

.. ipython:: python
    :suppress:

    import numpy as np
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    np.set_printoptions(precision=4, suppress=True)

single action with single trajectory
------------------------------------
.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)

multiple actions with multiple trajectories
-------------------------------------------
.. ipython:: python
    
    from pytraj.parallel import PJob

    tasklist = []
    tasklist.append((pt.radgyr, traj))
    tasklist.append((pt.molsurf, traj, '@CA'))

    # perform each action on each CPUs (total 2 CPUs)
    pjob = PJob(tasklist)
    print(pjob.compute())
