Installation
============

.. contents::

Requires
--------
if you install ``pytraj`` via ``conda``, conda will automatically install those libraries for you.

`numpy <http://www.numpy.org/>`_

`cython <http://cython.org/>`_, version >= 0.23

``python 2.7, 3.4, 3.5``

Install
-------

from conda
~~~~~~~~~~

We higly recommend install ``pytraj`` by `conda <http://conda.pydata.org/docs/intro.html>`_

.. code-block:: bash

    conda install -c ambermd pytraj-dev libcpptraj-dev

This takes only less than 30 seconds.

If you don't want to read long description about installing conda, just copy and paste below script to your terminal (for Linux).
For Mac user, you need to follow ``conda`` website.

.. code-block:: bash

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ sh Miniconda3-latest-Linux-x86_64.sh


from source code
~~~~~~~~~~~~~~~~

Alternatively, user can install ``pytraj`` from source code in `github``::

    git clone https://github.com/Amber-MD/pytraj
    cd pytraj
    python ./setup.py install

Depend on your machine, the fresh installation (``libcpptraj`` + ``pytraj``) could take 2 to 4 minutes.

Update pytraj
-------------

from conda
~~~~~~~~~~
If you install ``pytraj`` by conda, you can update it easily

.. code-block:: bash

    conda update -c ambermd pytraj-dev libcpptraj-dev

from github 
~~~~~~~~~~~

if you install from source code and want to update the development code in github, try to
follow below.

.. code-block:: bash
    
    $ # make sure to go to pytraj folder (which has README.md, ./tests ...)
    $ git pull
    $ python ./setup.py install

if you install ``pytraj`` via github and want to update ``cpptraj``

.. code-block:: bash

    $ cd cpptraj
    $ git pull
    $ make libcpptraj

Uninstall
---------

from conda
~~~~~~~~~~

.. code-block:: bash

    $ conda remove pytraj-dev libcpptraj-dev
    

Troubleshooting
---------------

``python setup.py install`` is standard process to install a new Python package.
But if you are new to Python and got ``permission denied`` error, try to install ``pytraj`` in your home folder.

.. code-block:: bash
    
    # install pytraj in $HOME/.local
    python ./setup.py install --user

    # or install pytraj in ``your_favorite_dir``
    python ./setup.py install --prefix=your_favorite_dir
    # if you do this, make sure to add ``your_favorite_dir`` to $PYTHONPATH 
    export PYTHONPATH=your_favorite_dir:$PYTHONPATH

If you want to see further options, check ``python setup.py install --help``

install anaconda for all python packages
----------------------------------------

we highly recommend to install anaconda that has all good python packages (``cython``, ``numpy``, ``sklearn``, ``pandas``, ...). Anaconda is totally free.

+ go to `its website <http://continuum.io/downloads#py34>`_, choose your platform and
  python version. It's better to pick up Python3
+ download file: For example, we downloaded ``Anaconda3-2.1.0-Linux-x86_64.sh`` (Python3
  version)
+ just run ``bash Anaconda3-2.1.0-Linux-x86_64.sh`` and follow instruction. That's it, you have a Python eco-system here.
