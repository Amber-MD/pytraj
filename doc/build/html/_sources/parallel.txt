.. _parallel:

Parallel calculation
====================

| 

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

Supported

.. contents::

.. ipython:: python
    :suppress:

    import numpy as np
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    np.set_printoptions(precision=4, suppress=True)

Disclaimer
----------

This is experimental design, syntax might be changed.

single action with single trajectory
------------------------------------

``pytraj`` use `python multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_, so users don't need to install extra package.

Example: parallel calculation with single action
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    pt.pmap(pt.radgyr, traj, n_cores=4)
    # auto-join the data
    pt.pmap(pt.radgyr, traj, n_cores=4, dtype='dict')

Example: parallel calculation with cpptraj's command style
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
    :suppress:

    import pytraj as pt
    from pytraj import matrix, vector, nmr
    from itertools import chain
    method_list_pmap = []
    method_list_openmp = []

    for method_str in chain(dir(pt), dir(matrix), dir(vector), dir(nmr)):
        try:
            method = getattr(pt, method_str)
            if hasattr(method, '_is_parallelizable') and method._is_parallelizable:
                method_list_pmap.append(method)
            if hasattr(method, '_openmp_capability') and method._openmp_capability:
                method_list_openmp.append(method)
        except AttributeError:
            pass

    pmap_ = []
    for method in set(method_list_pmap):
        name = str(method).split()[1]
        #if 'calc_' in name:
        #    name = name.split('calc_')[-1]
        pmap_.append(name)
    supported_pmap_methods = sorted(pmap_)

    openmp_ = []
    for method in set(method_list_openmp):
        name = str(method).split()[1]
        #if 'calc_' in name:
        #    name = name.split('calc_')[-1]
        openmp_.append(name)
    supported_openmp_methods = sorted(openmp_)


.. ipython:: python

    for method in supported_pmap_methods:
        print(method)

Supported methods for ``pmap`` if using cpptraj's command style
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**coming soon**


Supported methods for ``openmp``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

    for method in supported_openmp_methods:
        print(method)
    print("")


Rule of thumb for choosing ``pmap`` or ``openmp``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Always try to install ``pytraj`` and ``cpptraj`` with ``-openmp`` flag.
If method supports openmp, use openmp.

pmap doc
~~~~~~~~

:ref:`API <pytraj.pmap>`
