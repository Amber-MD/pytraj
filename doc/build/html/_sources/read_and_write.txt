.. _read_and_write:

Read and write
==============

(coming soon)

.. ipython:: python
    
    # read amber format
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    traj

    # write to CHARMM dcd format
    traj.save('test.dcd', overwrite=True)

    # write more than two trajs
    pt.write_traj('test2.dcd', [traj[:3], traj[5:]], top=traj.top, overwrite=True)

:ref:`trajectoryiterator` is able to load a list of files, even the total size does not fit into 
memory

.. code-block:: python

    In [1]: ls remd.x.*
    remd.x.000  remd.x.003  remd.x.006  remd.x.009  remd.x.012  remd.x.015
    remd.x.001  remd.x.004  remd.x.007  remd.x.010  remd.x.013  remd.x.016
    remd.x.002  remd.x.005  remd.x.008  remd.x.011  remd.x.014  remd.x.017
    
    In [2]: traj = pt.iterload('remd.x.*', 'myparm.parm7')
    
    In [3]: traj
    Out[3]:
    <pytraj.TrajectoryIterator, 40000 frames, include:
    <Topology: 17443 atoms, 5666 residues, 5634 mols, PBC with box type = truncoct>>
    
    In [4]: traj._estimated_MB
    Out[4]: 15969.54345703125
