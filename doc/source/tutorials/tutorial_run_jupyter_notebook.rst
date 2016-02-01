.. _tutorial_run_jupyter_notebook:

How to run Jupyter notebook
---------------------------

- Install::

   # if you already install python via Anaconda distribution
   # https://www.continuum.io/downloads, you don't need to install jupyter notebook again.
   conda install jupyter notebook

- Download notebook::
 
   wget  http://amber-md.github.io/pytraj/latest/tutorials/pca_evaluated.ipynb 

- Run::

   # tz2.{nc, parm7} can be downloaded from
   # https://github.com/Amber-MD/pytraj/tree/master/tests/data

   ipython notebook pca_evaluated.ipynb

   # or use AMBER16 distribution
   amber.ipython notebook pca_evaluated.ipython

- See also::

   Run Jupyter notebook `remotely <remote_jupyter_notebook>`_

- Example:

.. image:: ../images/tutorial_autoimage.png
