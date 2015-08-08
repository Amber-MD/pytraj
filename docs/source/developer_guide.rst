Developer guide for pytraj
==========================
(draft)

Python style guide
------------------
Try to follow `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_

Try to read `numpy doc <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

Python 2 and 3 compat
---------------------
Use `six <http://pythonhosted.org/six/>`_ to write your compat code. 

Add new method to pytraj
------------------------
Check ``pytraj.common_actions`` for example.

Test your code
--------------
New method, new change must have testing code.

Currently, all testing codes are in ``pytraj/tests/`` folder. We can use ``cp template_unittest.py test_new_method_name_example.py``. To run all tests, ``python ./run_all_and_find_fails.py``. To run tests having specific keywords ``python ./run_tests_with_keyword.py your_key_word``. Outputs from test scripts are saved to ``output.txt`` and error status is saved to ``log`` file.

The script ``./run_all_and_find_fails.py`` only look for file starting with ``test_`` and having key word ``unittest``. Check ``tests/get_unittest_files.py`` for further detail.

We're really happy to accept PR to update test, using `nosetests <https://nose.readthedocs.org/en/latest/>`_, `pytest <http://pytest.org/latest/>`_ or whatever reasonable.

External codes
--------------
Try to put all external codes (``six.py``, ...) in ``pytraj/externals/`` folder.

Licence info
------------
``pytraj`` always welcomes code contribution. It's recommended to put your name in the code you write. However, for the sake of clearness, just put something very short, like ``Copyright (c) 2010-2013 your_first_and_last_name`` and give full details of your contribution, license in ``pytraj/licenses/`` folder.

cython
------
It's recommended to use ``cython`` to write or wrap high performance code. Please don't use ``cimport numpy``, use `memoryview <http://docs.cython.org/src/userguide/memoryviews.html>`_ instead

Since ``pytraj`` will be bundled with AmberTools in Amber, it's important that we should commit cythonized file too. The main idea is that user only need C++ compiler and ``cpptraj``, nothing else.

Read Also
---------
`cpptraj-dev guide <https://github.com/mojyt/cpptraj/blob/master/doc/CpptrajDevlopmentGuide.lyx>`_

`pandas contributing guide <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_
