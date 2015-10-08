Developer guide for pytraj
==========================

.. contents::

General philosophy
------------------

Try your best to follow instruction below. If you have questions, do not hesitate to ask. Don't be afraid that your
code does not look right or pretty. We can discuss by openning an `issue <https://github.com/Amber-MD/pytraj/issues>`_

If having any suggestions, open an issue too.

Our github repo
---------------

How to contribute code?

Please read general instruction about git control in `pandas website
<http://pandas.pydata.org/pandas-docs/stable/contributing.html#version-control-git-and-github>`_

Our git repo is here:

.. code-block:: bash

    $ git clone https://github.com/Amber-MD/pytraj/
    $ cd pytraj
    $ git branch a_new_feature_or_whatever_name
    $ git checkout a_new_feature_or_whatever_name
    $ # do any work on this branch and make pull request

Python style guide
------------------
Try to follow `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_

Try to read `numpy doc <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

I (Hai) often use `yapf <https://github.com/google/yapf>`_ to format my code.

Python 2 and 3 compat
---------------------
Use `six <http://pythonhosted.org/six/>`_ to write your compat code. 
We put all common stuff in `pytraj.compat <https://github.com/Amber-MD/pytraj/blob/master/pytraj/compat.py>`_

.. note:: currently, I (Hai) are only working on Python3.

Install
-------
To speed up installation, please try to build in parallel.

.. code-block:: bash

    $ python ./setup.py build faster
    $ python ./setup.py install
    # if not see faster building or use only 1 core, use Ctrl-C and retry

Add new method to pytraj
------------------------
Check ``pytraj.common_actions`` for example.

Write your code for pytraj's parallel pmap
------------------------------------------

.. code-block:: python
 
    def new_method(traj, ...):
        #  make sure to use frame iterator like below
        for frame in traj:
            do_something_cool_with_frame
        return something_you_want

    # that's it. Now you can plug your method to ``pytraj.pmap``
    from pytraj import pmap
    pmap(n_cores=4, func=new_method, traj=traj, ...)

if you don't want to write code for `pmap`, just tag it with `noparallel` decorator

.. code-block:: python
    
    from pytraj.decorators import noparallel

    @noparallel
    def new_method(...):
        ...

Test your code
--------------
New method, new change must have testing code.

Currently, all testing codes are in **pytraj/tests/** folder. 

.. code-block:: bash

    # We can use
    $ cp template_unittest.py test_new_method_name_example.py
    # To run all tests
    $ python ./run_all_and_find_fails.py
    # To run tests having specific keywords 
    $ python ./run_tests_with_keyword.py your_key_word

Outputs from test scripts are saved to **output.txt** and error status is saved to **log** file.

The script ``./run_all_and_find_fails.py`` only look for file starting with ``test_`` and having key word ``unittest``. Check ``tests/get_unittest_files.py`` for further detail.

We're really happy to accept PR to update test, using `nosetests <https://nose.readthedocs.org/en/latest/>`_, `pytest <http://pytest.org/latest/>`_ or whatever reasonable.

External codes
--------------
Try to put all external codes (``six.py``, ...) in ``pytraj/externals/`` folder.

Licence info
------------
``pytraj`` always welcomes code contribution. It's recommended to put your name in the code you write. However, for the sake of clearness, just put something very short, like ``Copyright (c) 2010-2013 your_first_and_last_name`` and give full details of your contribution, license in ``pytraj/licenses/`` folder.

Build doc
---------

.. note:: make sure to install pytraj, cpptraj, numpy, ipython, matplotlib, memory_profiler, psutil. Install `sphinx-bootstrap-theme too <https://github.com/ryan-roemer/sphinx-bootstrap-theme>`_

.. code-block:: bash
    
    $ git clone https://github.com/Amber-MD/pytraj
    $ cd pytraj
    $ git checkout gh-pages
    $ cd doc
    $ make html

There are some tricks:

- let ipython run your code in ``.rst`` file by adding ipython directive::

   .. ipython:: python

- let ipython run your notebook and automatically convert to html file, add notebook directive::

    .. notebook:: data/plot_rmsd_radgyr_correlation.ipynb
       :skip_exceptions:

- let's see other tricks in::

    source/tutorials/*rst

How to make a tutorial and include it in pytraj's website
---------------------------------------------------------

I (Hai) prefer to use ipython notebook to write tutorial and include it in website. sphinx will run the notebook, convert to html file, insert it in correct page. 
But let's start with different ways to make a tutorial. First, make sure to::

  $ git checkout gh-pages

- use ipython directive: you just write the code and sphinx will run it for you. check::

  $ source/tutorials/basic_examples.rst

- use ipython notebook directive: you just write the code and sphinx will run it for you. This approach will have more richful layout. check::

  $ doc/source/tutorials/plot_correlation_matrix.rst

- Two above approaches are performed on the fly when you make the doc. If you don't want to rerun your notebook, you can run once, convert it to html file and include it in rst file::

  $ ipython nbconvert --to html your_notebook_name.ipynb
  $ # check doc/source/tutorials/lysozyme_order_parameter_.rst
  $ # (basically you just need to use .. raw:: html directive)

Push pytraj and libcpptraj to anaconda.org after successful build on travis
---------------------------------------------------------------------------

.. note:: This 'push' is for those who have permision to log in to ambermd account on anaconda.org

- website: `anaconda.org/ambermd <https://anaconda.org/ambermd>`_

- install ``ruby`` (google how)

- install ``travis``::

  $ gem install travis

- install anaconda-client::

  $ conda install anaconda-client
 
- In your terminal, log in to anaconda account::

  $ anaconda login
  $ # just enter your username and password

- generate anaconda token to give travis permision to push data in ambermd channel in anaconda.org::

  $ git clone https://github.com/Amber-MD/pytraj
  $ cd pytraj
  $ # generate token
  $ TOKEN=$(anaconda auth --create --name MyToken) 
  $ echo $TOKEN

- need to use ``travis`` to encrypt our token::

  $ travis encrypt TRAVIS_TO_ANACONDA=secretvalue

- make code change, commit, push to github so travis can build pytraj and libcpptraj::

  $ # after successful build, travis will push to anaconda.org by below command
  $ anaconda -t $TRAVIS_TO_ANACONDA upload --force -u ambermd -p pytraj-dev $HOME/miniconda/conda-bld/linux-64/pytraj-dev-*
  $ # check devtools/travis-ci/upload.sh and .travis.yml files for implementation.

cython
------
We recommended to use ``cython`` to write or wrap high performance code. Please don't use ``cimport numpy``, use `memoryview <http://docs.cython.org/src/userguide/memoryviews.html>`_ instead
Since ``pytraj`` will be bundled with AmberTools in Amber, it's important that we should commit cythonized file too. The main idea is that user only need C++ compiler and ``cpptraj``, nothing else.

For some unknow reasons, I (Hai) got segmentation fault if import numpy in the top of the module  when working with ``*.pyx`` file. It's better to import numpy locally (inside each method).


Read Also
---------
`cpptraj developer guide <https://github.com/mojyt/cpptraj/blob/master/doc/CpptrajDevlopmentGuide.pdf>`_

`test cpptraj api change with pytraj on travis <test_cpptraj_api>`_

`sklearn developer guide <http://scikit-learn.org/stable/developers/>`_

`pandas developer guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_
