.. _five_minutes:

Five minute learning
====================

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj


.. ipython:: python

    # load pytraj
    import pytraj as pt
    # load amber netcdf file
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj
    # perform analysis
    data = pt.molsurf(traj, mask='!@H=')
    data
    # plotting
    # write to disk in DCD format
    pt.write_traj('charmming.dcd', traj, overwrite=True)
    from matplotlib import pyplot as plt
    @savefig molsurf_plot.png
    plt.plot(data, '--bo')
