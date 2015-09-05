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

.. note:: make sure to install pytraj, cpptraj, numpy, ipython

.. code-block:: bash
    
    $ git clone https://github.com/Amber-MD/pytraj
    $ cd pytraj
    $ git checkout gh-pages
    $ cd doc
    $ make html

cython
------
It's recommended to use ``cython`` to write or wrap high performance code. Please don't use ``cimport numpy``, use `memoryview <http://docs.cython.org/src/userguide/memoryviews.html>`_ instead

Since ``pytraj`` will be bundled with AmberTools in Amber, it's important that we should commit cythonized file too. The main idea is that user only need C++ compiler and ``cpptraj``, nothing else.


Read Also
---------
`cpptraj-dev guide <https://github.com/mojyt/cpptraj/blob/master/doc/CpptrajDevlopmentGuide.pdf>`_

`Test cpptraj api change with pytraj on travis <test_cpptraj_api>`_

`pandas contributing guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_
