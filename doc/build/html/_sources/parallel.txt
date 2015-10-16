.. _parallel:

Parallel calculation
====================

.. contents::

.. note:: This is experimental design.

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj


``pytraj`` support different types of parallel computing.

``pytraj`` use `python multiprocessing <https://docs.python.org/3/library/multiprocessing.html>`_, so users don't need to install extra package.

.. ipython:: python
    :suppress:

    import numpy as np
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    np.set_printoptions(precision=4, suppress=True)

single action with single trajectory
------------------------------------

Example
~~~~~~~

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)

Supported methods for ``pmap``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
    :suppress:

    import pytraj as pt
    from pytraj import matrix, vector
    from itertools import chain
    method_list_pmap = []
    method_list_openmp = []

    for method_str in chain(dir(pt), dir(matrix), dir(vector)):
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
        if 'calc_' in name:
            name = name.split('calc_')[-1]
        pmap_.append(name)
    supported_pmap_methods = sorted(pmap_)

    openmp_ = []
    for method in set(method_list_openmp):
        name = str(method).split()[1]
        if 'calc_' in name:
            name = name.split('calc_')[-1]
        openmp_.append(name)
    supported_openmp_methods = sorted(openmp_)


.. ipython:: python

    for method in supported_pmap_methods:
        print(method)


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

multiple actions with multiple trajectories
-------------------------------------------

Only works with Python 3.

.. ipython:: python
    
    from pytraj.parallel import PJob

    tasklist = []
    tasklist.append((pt.radgyr, traj))
    tasklist.append((pt.molsurf, traj, '@CA'))

    # perform each action on each CPUs (total 2 CPUs)
    pjob = PJob(tasklist)
    print(pjob.compute())
