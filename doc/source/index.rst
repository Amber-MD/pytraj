.. pytraj documentation master file, created by
   sphinx-quickstart on Mon Jun 22 19:09:20 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pytraj's documentation - nothing to see yet!
=======================================================

**Overview**

``pytraj`` is a Python front-end of the popular ``cpptraj`` package. Its aim is to expose
``cpptraj``'s funtions to Python's ecosystem. Enjoy.

**Contents**

.. toctree::
   :maxdepth: 1

   installation
   basic_examples
   trajectory
   analysis
   modify_traj
   atom_mask_selection
   parallel
   whatsnew
   faq
   developer_guide
   misc
   cookbook
   api


**Simple code**

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    print(traj)
    data = pt.rmsd(traj, ref=0, mask='@CA')
    print(data)

**Contact**

* For bug reports, feature requests, code discussion ... please use the `GitHub issue tracker <https://github.com/Amber-MD/pytraj/issues>`_
* For scientific discussion please send email to `amber <http://ambermd.org/#correspondence>`_
* For things that directly related to ``cpptraj``, please open `issue here
  <https://github.com/Amber-MD/cpptraj/issues>`_

**Indices and tables**

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
