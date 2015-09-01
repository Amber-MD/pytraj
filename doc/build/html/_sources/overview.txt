.. _overview:

Overview
========

**pytraj** is a Python front-end package of the popular `cpptraj <http://pubs.acs.org/doi/abs/10.1021/ct400341p>`_ program.

#. Why bother using **pytraj**? 

    If you want to extend cpptraj's functionanity to Python eco-system such as `numpy <http://www.numpy.org/>`_, `scipy <http://www.scipy.org/>`_, `pandas <http://pandas.pydata.org/>`_, `sklearn <http://scikit-learn.org/stable/>`_ and many more. Check some examples `here <tutorials/mdtraj_adapted>`_

#. What does **pytraj** offer? 

   + `support many types of analysis <analysis>`_
   + `support parallel processing <parallel>`_
   + `able to handle many files at the same time <process_many_files>`_
   + `able to handle very large trajectory <design_trajectory>`_
   + simple (and fast) calculation

      import pytraj as pt
      traj = pt.iterload('traj.nc', 'my_parm.top')
      pt.radgyr(traj)

   + include many cpptraj features
       + automatic file detection by its **content**
       + smart ``autoimage`` function for PBC simulation.
       + ... (to be filled soon) 
   + `able to be used in other langues, such as Julia <julia>`_

#. How to get started

   + `check installation <installation>`_
   + do some :ref:`tutorials`
   + get reference by :ref:`genindex`
