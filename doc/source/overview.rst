.. _overview:

Overview
========

**pytraj** is a Python package as a front-end package of the popular `cpptraj <http://pubs.acs.org/doi/abs/10.1021/ct400341p>`_ program.

#. Why bother using **pytraj**? If you want to extend cpptraj's functionanity to Python eco-system such as `numpy <http://www.numpy.org/>`_, `scipy <http://www.scipy.org/>`_, `pandas <http://pandas.pydata.org/>`_, `sklearn <http://scikit-learn.org/stable/>`_ and many more

#. What does **pytraj** offer? 
   `support many types of analysis <analysis>`_
   `support parallel processing <parallel>`_
   `able to handle many files at the same time <read_and_write>`_
   `able to handle very large trajectory <design_trajectory>`_
   simple (hopefully) usage ::

      import pytraj as pt
      traj = pt.iterload('traj.nc', 'my_parm.top')
      pt.radgyr(traj)

#. How to get started
   `check installation <installation>`_
   do some tutorials :ref:`tutorial`
