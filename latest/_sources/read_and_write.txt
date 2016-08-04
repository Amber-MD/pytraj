.. _read_and_write:


Supported file formats
======================

.. contents::

Read
----

Most commmon usage: load a single file with all frame

.. ipython:: python
    
    # read amber format
    import pytraj as pt
    traj = pt.iterload('tz2.ortho.nc', 'tz2.ortho.parm7')
    traj

Most commmon usage: load a single file with frame stride

.. code-block:: python
    
    # load only frame 0, 2, 4, 6 (python use 0-based index and skip final index)
    traj = pt.iterload('tz2.nc', 'tz2.parm7', frame_slice=[(0, 8, 2),])

Load many files at once

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

Load many files with frame stride

Example: load 1000 frames from two trajectories (500 each), skip every two frames.

.. code-block:: python

    In [1]: traj = pt.iterload(['remd.x.000', 'remd.x.001'], 'myparm.parm7', frame_slice=[(0, 1000, 2), (0, 1000, 2)])
    
    In [2]: traj
    Out[2]: 
    <pytraj.TrajectoryIterator, 1000 frames, include:
    <Topology: 17443 atoms, 5666 residues, 5634 mols, PBC with box type = truncoct>>

which is similiar to ``cpptraj`` input:

.. code-block:: bash

    parm myparm.parm7
    trajin remd.x.000 1 1000 2
    trajin remd.x.001 1 1000 2

.. note:: cpptraj uses 1-based index for input while ``pytraj`` used 0-based index.


Write
-----

.. ipython:: python

    # write to CHARMM dcd format
    pt.write_traj('test.dcd', traj, overwrite=True)

    # write with given frames
    pt.write_traj('test2.dcd', traj, frame_indices=[0, 3, 7, 9], overwrite=True)

Supported file formats
----------------------

.. autofunction:: pytraj.write_traj
