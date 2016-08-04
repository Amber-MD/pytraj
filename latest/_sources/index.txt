.. pytraj documentation master file, created by
   sphinx-quickstart on Mon Jun 22 19:09:20 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

.. raw:: html
    :file: _static/track.html

Welcome
=======

``pytraj`` is a Python front-end of the popular ``cpptraj`` package. Its aim is to expose
``cpptraj``'s funtions to Python's ecosystem. Enjoy.

.. raw:: html

   <div class="col-md-3">
   <h2>Overview</h2>

.. toctree::
   :maxdepth: 2

   overview

`AMBER16 users - click me <http://amber-md.github.io/pytraj/release/1.0.4/index.html>`_

`Release versions <http://amber-md.github.io/pytraj/release/>`_

.. raw:: html
    :file: versions.html

.. raw:: html

   </div>
   <div class="col-md-3">
   <h2>Documentation</h2>

.. toctree::
   :maxdepth: 1

   installation
   trajectory
   topology
   read_and_write
   tutorials/index
   analysis
   modify_traj
   atom_mask_selection
   trajectory_slice
   parallel
   whatsnew
   faq
   developer_guide
   misc
   conda
   cookbook
   trajectory_viewer
   api

.. raw:: html

   </div>
   <div class="col-md-3">
   <h2>Plot</h2>

.. image:: images/PCA_heart.png
   :target: tutorials/plot.html

|

**Indices and tables**

* `fork and contribute <https://github.com/Amber-MD/pytraj>`_
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. raw:: html

   </div>
   </div>
   </div>

.. raw:: html

   </div>
   <div class="col-md-3">
   <h3>Jupyter notebook</h3>

.. image:: http://jupyter.org/assets/jupyterpreview.png
   :target: http://jupyter.org/
   :height: 200


.. raw:: html

   </div>
   <div class="col-md-3">
   <h3><a href=trajectory_viewer.html> Trajectory visualization </a></h3>

.. raw:: html
   :file: ngl_example.html

.. raw:: html

   </div>
   <div class="col-md-3">
   <h3>Try pytraj online</h3>

.. image:: http://mybinder.org/assets/images/logo.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. raw:: html

   </div>
