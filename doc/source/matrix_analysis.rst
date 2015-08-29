.. _matrix_analysis:

Matrix Analysis
===============

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

Examples
--------

.. code-block:: python
    
    from pytraj import matrix_analysis as ma
    ma.distance_matrix(traj, '@CA')
    ma.idea_matrix(traj, '@CA')
    ma.correlation_matrix(traj, '@CA')

doc
---

.. automodule:: pytraj.matrix_analysis
