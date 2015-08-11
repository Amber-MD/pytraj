Parallel calculation
--------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)
