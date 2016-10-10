.. _parallel:

Parallel calculation
====================

| 

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

Supported

.. note:: pytraj does not support writing trajectory in parallel. Please check cpptraj manual for this option.

.. contents::

.. ipython:: python
    :suppress:

    import numpy as np
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    np.set_printoptions(precision=4, suppress=True)

.. note:: This is experimental design, syntax might be changed.


Example: parallel calculation with single action
------------------------------------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    pt.pmap(pt.radgyr, traj, n_cores=4)
    # auto-join the data
    pt.pmap(pt.radgyr, traj, n_cores=4)

Example: parallel calculation with cpptraj's command style
----------------------------------------------------------

**Rule of thumb**: Only support limited cpptraj's actions (not analysis). For now, there
is no list of supported actions yet but you can perform any action in parallel if the
result does not depend on the number of frame. For example, you can NOT use 'atomicfluct'
or 'matrix covar', ... This parallel method is GOOD for distance, angle, rms, dihedral,
vector, multidihedral, ...

.. ipython:: python

    # perform autoimage, then rmsfit to 1st frame
    # create a reference
    ref = traj[0]
    ref
    # need to explicitly do autoimage to reproduce cpptraj's serial run
    # since `ref` is a single Frame, we need to provide Topology for `pt.autoimage`
    pt.autoimage(ref, top=traj.top)
    # add to reference list
    reflist = [ref, ]
    # need to explicitly specify reference by using 'refindex'
    # 'refindex 0' is equal to 'reflist[0]' and so on.
    # we need to explicitly specify the index so pytraj can correctly send the reference to different cores.
    pt.pmap(['autoimage', 'rms refindex 0', 'distance :3 :8', 'multidihedral phi psi resrange 3-5'], traj, ref=reflist, n_cores=4)

    # compare to cpptraj's serial run
    state = pt.load_cpptraj_state('''
    parm tz2.ortho.parm7
    trajin tz2.ortho.nc
    autoimage
    rms
    distance :3 :8
    multidihedral phi psi resrange 3-5
    ''')
    state.run()
    # pull output from state.data
    # since cpptraj store Topology in `state.data`, we need to skip it
    state.data[1:].to_dict()


Supported methods for ``pmap`` if using pytraj's methods
--------------------------------------------------------

.. ipython:: python

    import pytraj as pt
    for method in pt.utils.misc.parallel_info('pmap'):
        print(method)

Supported methods for ``pmap`` if using cpptraj's command style
---------------------------------------------------------------

**coming soon**


Supported methods for ``openmp``
--------------------------------

.. ipython:: python

    for method in pt.utils.misc.parallel_info('openmp'):
        print(method)
    print("")


Rule of thumb for choosing ``pmap`` or ``openmp``?
--------------------------------------------------

Always try to install ``pytraj`` and ``cpptraj`` with ``-openmp`` flag.
If method supports openmp, use openmp.

Benchmark
---------

System info::

    format: AMBER netcdf file

    pytraj.TrajectoryIterator, 200000 frames: 
    Size: 58.150291 (GB)
    <Topology: 13008 atoms, 4189 residues, 4174 mols, PBC with box type = truncoct>

    method: pytraj.rmsd (please check the script below)


Multiprocessing
~~~~~~~~~~~~~~~

.. image:: images/bench_pmap_casegroup.png


Distributed: using MPI (via mpi4py)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: images/bench_pmap_mpi.png

Script for multiprocessing

.. code-block:: python

    from multiprocessing import cpu_count
    from glob import glob
    import pytraj as pt
    from time import time
    
    print("max cores = {}".format(cpu_count()))
    filenames = glob('mdx/md*.nc') * 10
    traj = pt.iterload(filenames, 'prmtop')
    print(traj)
    
    mask = '!:WAT,Na+'
    
    t0 = time()
    pt.rmsd(traj, mask=mask)
    t_serial = time() - t0
    # print('serial time = ', t_serial)
    
    func = pt.rmsd
    
    # for n_cores in [1, 2, 4, 6, 8, 16, 20, 24]:
    for n_cores in [1, 2, 4, 6, 7]:
        t0 = time()
        pt.pmap(func, traj, mask=mask, n_cores=n_cores, ref=traj[0])
        t_par = time() - t0
        print(n_cores, t_serial / t_par)

Script for MPI

.. code-block:: python

    fromm glob import glob
    import pytraj as pt
    from time import time
    from mpi4py import MPI
    
    comm = MPI.COMM_WORLD
    
    # need to run program in serial and update serial time
    serial_time = 45.
    
    filenames = glob('mdx/md*.nc') * 10
    traj = pt.iterload(filenames, 'prmtop')
    
    mask = '!:WAT,Na+'
    
    func = pt.rmsd
    
    t0 = time()
    x = pt.pmap_mpi(func, traj, mask=mask, ref=traj[0])
    t_par = time() - t0
    
    if comm.rank == 0:
        print(serial_time/t_par)


pmap API doc
------------

:ref:`API <pytraj.pmap>`
