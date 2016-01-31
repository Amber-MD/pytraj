.. _tips:

Tips
====

.. contents::

.. _process_many_files:

process many files
------------------

Normally, user needs to follow below code to process many files

.. ipython:: python

    import pytraj as pt
    import numpy as np
    template = 'tz2.%s.nc'
    flist = [template % str(i) for i in range(4)]
    print(flist)
    trajlist = [pt.iterload(fname, 'tz2.parm7') for fname in flist]

    # calculate radgyr
    data = []
    for traj in trajlist:
        data.append(pt.radgyr(traj))
    data = np.array(data).flatten()
    print(data)

However, ``pytraj`` offer shorter (and easy?) way to do

.. ipython:: python
    
    # load all files into a single TrajectoryIterator
    filelist = ['tz2.0.nc', 'tz2.1.nc', 'tz2.2.nc']
    traj = pt.iterload(filelist, 'tz2.parm7')
    # perform calculation
    data = pt.radgyr(traj)
    # that's it
    print(data)

``pytraj`` will auto-allocate ``data`` too

.. ipython:: python
    
    print([t.filename.split('/')[-1] for t in trajlist])
    data = pt.radgyr(trajlist, top=trajlist[0].top)
    print(data)

load from a list of files with frame stride
-------------------------------------------

Supposed you have a list of 5 (or whatever) trajectories, you only want to load those files 1st to 100-th frames
and skip every 10 frames. Below is a convention ``cpptraj`` input.

.. code-block:: bash

    parm 2koc.parm7
    trajin traj0.nc 1 100 10
    trajin traj1.nc 1 100 10
    trajin traj2.nc 1 100 10
    trajin traj3.nc 1 100 10
    trajin traj4.nc 1 100 10

In ``pytraj``, you can specify ``frame_slice``

.. code-block:: python

    import pytraj as pt
    pt.iterload('traj*.nc', top='2koc.parm7', frame_slice=[(0, 100, 10),]*5)

    # [(0, 100, 10),]*5 is equal to [(0, 100, 10), (0, 100, 10),(0, 100, 10),(0, 100, 10),(0, 100, 10),]

load specific frame numbers to memory
-------------------------------------

.. ipython:: python

    import pytraj as pt
    frame_indices = [2, 4, 7, 51, 53]
    # use ``load`` to load those frames to memory
    traj0 = pt.load('tz2.nc', 'tz2.parm7', frame_indices=frame_indices)
    traj0

    # only loadd coordinates for specific atoms
    traj1 = pt.load('tz2.nc', 'tz2.parm7', frame_indices=frame_indices, mask='@CA')
    traj1

    # or use ``iterload``
    frame_indices = [2, 4, 7, 51, 53]
    traj2 = pt.iterload('tz2.nc', 'tz2.parm7')
    traj2
    traj2[frame_indices, '@CA']

merge multiple trajectories to a single file
--------------------------------------------

.. ipython:: python

    import pytraj as pt
    # load multiple files
    traj = pt.iterload(['tz2.0.nc', 'tz2.1.nc', 'tz2.2.nc'], top='tz2.parm7')
    traj.save('tz2_combined.nc', overwrite=True)

memory saving
-------------

If memory is critical, do not load all frames into memory.

.. ipython:: python

    # DO this (only a single frame will be loaded to memory)
    pt.radgyr(traj, frame_indices=[0, 200, 300, 301])

    # DON'T do this if you want to save memory (all 4 frames will be loaded to memory)
    pt.radgyr(traj[[0, 200, 300, 301]])

    pt.iterframe(traj, frame_indices=[0, 200, 300, 301])
    traj[[0, 200, 300, 301]]


Example: calculate pairwise rmsd for 3000 frames (~3.9 G) only costs 112.7 MB

.. code-block:: bash

    $ python memory/pairwise_rmsd_shorter_version.py

    pytraj.TrajectoryIterator, 3000 frames:
    Size: 3.911264 (GB)
    <Topology: 58329 atoms, 19202 residues, 19168 mols, PBC with box type = ortho>
    
    Filename: memory/pairwise_rmsd_shorter_version.py
    
    Line #    Mem usage    Increment   Line Contents
    ================================================
         6     30.1 MiB      0.0 MiB   @profile
         7                             def compute(mask='!:WAT,K+,Cl-', mat_type='half'):
         8     64.0 MiB     33.9 MiB       traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo', frame_slice=(0, 3000))
         9     64.1 MiB      0.1 MiB       print(traj)
        10    112.7 MiB     48.5 MiB       return pt.pairwise_rmsd(traj, mask=mask, mat_type=mat_type)


See also: :ref:`trajectory_slice`

convert trajectory
------------------

.. code-block:: python
    
    # convert Amber netcdf to Charmm dcd file.
    pt.iterload('traj.nc', 'prmtop').save('traj.dcd', overwrite=True)

pickle data
-----------

Sometimes you need to perform very long analysis (hours), you need to save the output to
disk to do further analysis. You have options to save data to different files and write
code to load the data back. However, you can use ``pytraj.to_pickle`` and
``pytraj.read_pickle`` to save the state of data. Check the example:

.. ipython:: python

    traj3 = pt.load_pdb_rcsb('1l2y')
    data = pt.dssp(traj, ':3-7')
    data
    pt.to_pickle(data, 'my_data.pk')
    # load the data's state back for further analysis
    pt.read_pickle('my_data.pk')
    # note: do not read_pickle from files that don't belong to you. It's not secure.

speed up calculation with parallel: using MPI
---------------------------------------------

.. code-block:: bash

    $ cat > radgyr_mpi.py <<EOF
    import pytraj as pt
    
    # import mpi4py to get rank of each core
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    
    # load trajectory to each core. Use iterload to save memory
    traj = pt.iterload('tz2.nc', 'tz2.parm7')

    # compute radgyr by sending this method to pt.pmap_mpi function
    data = pt.pmap_mpi(pt.radgyr, traj, '@CA')
    
    # data is sent to first core (rank=0)
    if comm.rank == 0:
        # save data
        # pt.to_pickle(data, 'data.pk')
        print(data)
    EOF

    $ # run
    $ mpirun -n 4 python radgyr_mpi.py
    # you can reload data by `pt.read_pickle('data.pk')

speed up calculation with parallel: multiple cores
--------------------------------------------------

.. ipython:: python
    
    # send pt.radgyr method to `pt.pmap` function
    # need to specify the number of cores
    # (choose wisely)

    pt.pmap(pt.radgyr, traj, '@CA', n_cores=4)

 
read cpptraj manual
-------------------

This does not work with ipython-notebook but it's still good for interactive ipython

.. code-block:: python

    In [106]: import pytraj as pt
    In [107]: pt.info('radgyr')
            [<name>] [<mask1>] [out <filename>] [mass] [nomax] [tensor]
              Calculate radius of gyration of atoms in <mask>
