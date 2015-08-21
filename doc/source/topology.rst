.. _topology:

Topology
========

Load topology
-------------

From **iterload**

.. code-block:: python

    import pytraj as pt
    traj = pt.iterload('traj.nc', '2KOC.parm7')
    top = traj.top

From **load_topology**

.. code-block:: python

    top = pt.load_topology('2KOC.parm7')

Supported file format
---------------------

.. todo:: add table

Topology editing
----------------
pytraj and cpptraj offer limited Topology editing. For more features and power, please do
visit `ParmEd <http://parmed.github.io/ParmEd/html/index.html>`_

* slicing a Topology
* save to different format
