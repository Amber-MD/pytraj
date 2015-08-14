.. _tips:

Tips
====

* turn on cpptraj's verbose

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
    traj
    pt.set_cpptraj_verbose(True)
    pt.radgyr(traj, '@C,N,O,H')
