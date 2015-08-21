.. _topology:

Topology
========

.. note:: pytraj and cpptraj offer limited Topology editing. For more features and power, please do
visit `ParmEd <http://parmed.github.io/ParmEd/html/index.html>`_

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
(table was adapted from Amber15 manual)

Format     Extension Notes
======     ========= =====
Amber      parm7     Write/Read
PDB        pdb       Read Only
Mol2       mol2      Read Only
CIF        cif       Read Only
Charmm PSF psf       Limited PSF Write
SDF        sdf       Read Only
Tinker ARC arc       Read Only

.. todo:: add table

Topology editing
----------------

* slicing a Topology
* save to different format
