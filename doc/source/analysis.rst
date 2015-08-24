.. _analysis:

Analysis
========

.. ipython:: python
    :suppress:

    import numpy as np
    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')

.. todo: make a script to write

.. seealso::

    `Methods that modify Trajectory coordinates <modify_traj>`_


- distance_
- angle_
- dihedral_ 
- rmsd_
- rmsd_perres_
- pairwise_rmsd_ 
- radgyr_ 
- molsurf_ 
- center_of_mass_
- center_of_geometry_
- dssp_ 
- multidihedral_ 
- bfactors_ 
- search_hbonds_
- autoimage_
- jcoupling_ 
- density_
- volume_
- minimum_distance_
- lifetime_
- auto_correlation_function_
- principal_axes_
- cross_correlation_function_
- timecorr_
- randomize_ions_
- crank_
- closest_
- search_neighbors_
- watershell_
- :ref:`nucleic_acid`


.. _distance:
.. autofunction:: pytraj.distance

.. _angle:
.. autofunction:: pytraj.angle

.. _dihedral:
.. autofunction:: pytraj.dihedral

.. _rmsd:
.. autofunction:: pytraj.rmsd

.. _rmsd_perres:
.. autofunction:: pytraj.rmsd_perres

.. _pairwise_rmsd:
.. autofunction:: pytraj.pairwise_rmsd

.. ipython:: python
    pt.pairwise_rmsd(traj, '@CA')

.. _dssp:
.. autofunction:: pytraj.dssp

.. _multidihedral:
.. autofunction:: pytraj.multidihedral

.. _bfactors:
.. autofunction:: pytraj.bfactors 

.. _radgyr:
.. autofunction:: pytraj.radgyr

.. _molsurf:
.. autofunction:: pytraj.molsurf

.. _center_of_mass:
.. autofunction:: pytraj.center_of_mass

.. _center_of_geometry:
.. autofunction:: pytraj.center_of_geometry

.. _search_hbonds:
.. autofunction:: pytraj.search_hbonds

.. _jcoupling:
.. autofunction:: pytraj.jcoupling

.. _density:
.. autofunction:: pytraj.density

.. _volume:
.. autofunction:: pytraj.volume

.. _minimum_distance:
.. autofunction:: pytraj.mindist

.. _lifetime:
.. autofunction:: pytraj.lifetime

.. _auto_correlation_function:
.. autofunction:: pytraj.acorr

.. _cross_correlation_function:
.. autofunction:: pytraj.xcorr

.. _timecorr:
.. autofunction:: pytraj.timecorr

.. _principal_axes:
.. autofunction:: pytraj.principal_axes

.. _randomize_ions:
.. autofunction:: pytraj.randomize_ions

.. _crank:
.. autofunction:: pytraj.crank

.. _autoimage:
.. autofunction:: pytraj.autoimage

.. _closest:
.. autofunction:: pytraj.closest

.. _search_neighbors:
.. autofunction:: pytraj.search_neighbors

.. _watershell:
.. autofunction:: pytraj.watershell
