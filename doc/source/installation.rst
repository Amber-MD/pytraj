Installation
============

.. include:: mybinder.rst

.. contents::

.. note:: For Linux user, we highly recommend install pytraj via `conda`

Supported platforms
-------------------
- Linux
- OSX

Supported Python versions
-------------------------
- 2.7, >=3.4

Requires
--------

- numpy

- (optional) `cython <http://cython.org/>`_, >= 0.23. Cython is only required if installing development version::

    pip install cython --upgrade

    # conda install cython



Install
-------

from conda (Linux, OX)
~~~~~~~~~~~~~~~~~~~~~~

We higly recommend install ``pytraj`` by `conda <http://conda.pydata.org/docs/intro.html>`_

.. code-block:: bash

    # install latest version
    conda install -c ambermd pytraj

    # install specific version
    conda install -c ambermd pytraj==1.0.8

This takes only less than 30 seconds.

If you don't want to read long description about installing conda, just copy and paste below script to your terminal (for Linux).
For Mac user, you need to follow ``conda`` website.

.. code-block:: bash

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ sh Miniconda3-latest-Linux-x86_64.sh

From pip (Linux, OSX)
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash
    
    # latest version
    pip install pytraj

    # specific version
    pip install pytraj==1.0.8


From source code (easy way: Linux, OSX)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/Amber-MD/pytraj
    cd pytraj

    # linux
    python ./setup.py install

    # osx:
    python setup.py install

    # note: pytraj will automatically install cpptraj first.

Depend on your machine, the fresh installation (``libcpptraj`` + ``pytraj``) could take 2 to 4 minutes.

From source code (hard way - expert only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you want to install `libcpptraj.so` by yourself.

- First, download cpptraj::

    git clone https://github.com/Amber-MD/cpptraj
    cd cpptraj
    bash configure -shared -openmp gnu

    # if you are AMBER user, you can add -amberlib
    bash configure -shared -openmp -amberlib gnu
    make libcpptraj -j4

    # please check bash configure --full-help for other options.
    # check: https://github.com/Amber-MD/cpptraj too
    export CPPTRAJHOME=`pwd`

- Then, install ``pytraj`` ::

    # cd to any folder you want to store pytraj code
    # then
    git clone https://github.com/Amber-MD/pytraj
    cd pytraj
    python ./setup.py install


From AMBER distribution (Linux, OSX)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pytraj is included in AMBER (version >= 16): ambermd.org


Update pytraj
-------------

Rule of thumb: using the same tool to install/update (upgrade)

From conda
~~~~~~~~~~
If you install ``pytraj`` by conda, you can update it easily

.. code-block:: bash

    conda update -c ambermd pytraj --force

From pip
~~~~~~~~

.. code-block:: bash

    pip install --upgrade pytraj


From github 
~~~~~~~~~~~

if you install from source code and want to update the development code in github, you
need to update both `libcpptraj` and `pytraj`

.. code-block:: bash

    $ cd /to/pytraj/root/folder
    $ cd cpptraj
    $ git pull
    $ make libcpptraj -j8

then

.. code-block:: bash
    
    $ cd /to/pytraj/root/folder
    $ git pull
    $ python ./setup.py install


Update to AMBERHOME
~~~~~~~~~~~~~~~~~~~


- From release version

.. code-block:: bash

    $ amber.pip install pytraj --prefix=$AMBERHOME
    # if you don't have amber.pip, just use pip

- From github development version

.. code-block:: bash

    rm $AMBERHOME/lib/libcpptraj.so # so we can use cpptraj github version 
    # (which will be installed by pytraj itself)

    cd /to/your/favorite/dir
    git clone https://github.com/Amber-MD/pytraj
    cd pytraj

    # if there's no git
    # wget https://github.com/Amber-MD/pytraj/archive/master.zip
    # unzip master
    # cd pytraj-master

    python setup.py install --prefix=$AMBERHOME # overwrite old version

    
Uninstall
---------

Rule of thumb: using the same tool to install/uninstall

From conda
~~~~~~~~~~

.. code-block:: bash

    $ conda remove pytraj libcpptraj

From pip
~~~~~~~~

.. code-block:: bash

    $ pip uninstall pytraj
    

Troubleshooting
---------------

Permission denied
~~~~~~~~~~~~~~~~~

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


Install ipython and its notebook for interactive data exploration
-----------------------------------------------------------------

`ipython <http://ipython.org/>`_ and its notebook is great program for interactive exloration of MD data.
Curious about how the notebook looks like? check out our `pairwise rmsd tutorial <http://amber-md.github.io/pytraj/doc/build/html/tutorials/tut_pairwise_rmsd.html>`_

If you are using ``anaconda``, just type ``ipython notebook``. If you have not haved ipython and its notebook, try ``conda install ipython``
For further instruction and information about ``ipython-notebook``, please check its website.

Install anaconda for all python packages
----------------------------------------

we highly recommend to install anaconda that has all good python packages (``cython``, ``numpy``, ``sklearn``, ``pandas``, ...). Anaconda is totally free.

+ go to `its website <http://continuum.io/downloads#py34>`_, choose your platform and
  python version. It's better to pick up Python3
+ download file: For example, we downloaded ``Anaconda3-2.1.0-Linux-x86_64.sh`` (Python3
  version)
+ just run ``bash Anaconda3-2.1.0-Linux-x86_64.sh`` and follow instruction. That's it, you have a Python eco-system here.


install and run jupyter noteook
-------------------------------

.. code-block:: bash

    # install
    conda install notebook

    # run
    jupyter notebook
    
    # or run
    jupyter notebook {your_notebook_name}.ipynb

If you want to run Jupyter notebook remotely, check :ref:`remote_jupyter_notebook`
