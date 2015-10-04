.. _overview:

Overview
========

.. contents::

**pytraj** is a Python front-end package of the popular `cpptraj <http://pubs.acs.org/doi/abs/10.1021/ct400341p>`_ program.

Why bother using **pytraj**? 
----------------------------

If you want to extend cpptraj's functionanity to Python eco-system such as `numpy <http://www.numpy.org/>`_, `scipy <http://www.scipy.org/>`_, `pandas <http://pandas.pydata.org/>`_, `sklearn <http://scikit-learn.org/stable/>`_ and many more. Check some examples `here <tutorials/mdtraj_adapted>`_

What does **pytraj** offer? 
---------------------------

+ `support many types of analysis <analysis>`_
+ `support parallel processing <parallel>`_
+ :ref:`able to handle many files at the same time <process_many_files>`
+ `able to handle very large trajectory <design_trajectory>`_
+ simple (and fast) calculation::

   import pytraj as pt
   traj = pt.iterload('traj.nc', 'my_parm.top')
   pt.radgyr(traj)

+ include many cpptraj features
    + automatic file detection by its **content**
    + smart ``autoimage`` function for PBC simulation.
    + ... (to be filled soon) 
+ `able to be used in other langues, such as Julia <julia>`_

How to get started
------------------

+ `check installation <installation>`_
+ `learn in 5 minutes <five_minutes>`_
+ `general concept <general_concept>`_
+ `check a list of things pytraj/cpptraj can do <analysis>`_

.. _contributors_and_citations:

Contributors
------------

`A list of contributors <https://github.com/Amber-MD/pytraj/blob/master/contributors/pytraj.rst>`_

Citations
---------

if you would like to acknowledge our works, please do consider citing both ``cpptraj`` and ``pytraj`` papers ::

    PTRAJ and CPPTRAJ : Software for Processing and Analysis of Molecular Dynamics Trajectory Data
    Daniel R. Roe and Thomas E. Cheatham, III
    Journal of Chemical Theory and Computation 2013 9 (7), 3084-3095 
    
::

    pytraj (https://github.com/Amber-MD/pytraj), Hai Nguyen,  et al. (2015) (in preperation)

.. note:: list of authors in ``pytraj`` is not complete. We will fill it soon.


Contact
-------

We always welcome your contribution, including (but not limited to) adding new code, bug report,
suggestion about website, doc, ... or anything that comes up. Feel free to contact us.
Even you do not like ``pytraj`` or some of its features, just tell us so we can make it
better. Cheers.

* For bug reports, feature requests, code discussion ... please use the `GitHub issue tracker <https://github.com/Amber-MD/pytraj/issues>`_
* For scientific discussion please send email to `amber <http://ambermd.org/#correspondence>`_
* For things that directly related to ``cpptraj``, please open `issue here
  <https://github.com/Amber-MD/cpptraj/issues>`_
