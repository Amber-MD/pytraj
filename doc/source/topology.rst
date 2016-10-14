.. _topology:

Topology
========

.. include:: mybinder.rst

.. note:: pytraj and cpptraj offer limited Topology editing. For more features and power, please do visit `ParmEd <http://parmed.github.io/ParmEd/html/index.html>`_

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

Supported file formats
----------------------

.. note:: cpptraj/pytraj (and parmed) recognizes file by its content, not by extention. So it's 'legal' to use ::

    pt.load_topology('2koc.very_very_long_ext')

(table was adapted from Amber15 manual)

========== ========= =================
Format     Extension Notes
========== ========= =================
Amber      parm7     Write/Read
PDB        pdb       Read Only
Mol2       mol2      Read Only
CIF        cif       Read Only
Charmm PSF psf       Limited PSF Write
SDF        sdf       Read Only
Tinker ARC arc       Read Only
========== ========= =================

Topology editing
----------------

* slicing a Topology ::

    top = pt.load_topology('2koc.pdb')
    top['@CA']

* save to different format ::

    top = pt.load_topology('2koc.parm7')
    top.save('test.mol2')
