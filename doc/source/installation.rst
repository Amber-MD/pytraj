Installation
============

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj


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

- python 2.7, >=3.4

- `cython <http://cython.org/>`_, >= 0.23. Cython is only required if installing development version::

    pip install cython --upgrade

    # conda install cython


Install
-------

from conda (Linux only)
~~~~~~~~~~~~~~~~~~~~~~~

We higly recommend install ``pytraj`` by `conda <http://conda.pydata.org/docs/intro.html>`_

.. code-block:: bash

    conda install -c ambermd pytraj-dev

This takes only less than 30 seconds.

If you don't want to read long description about installing conda, just copy and paste below script to your terminal (for Linux).
For Mac user, you need to follow ``conda`` website.

.. code-block:: bash

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ sh Miniconda3-latest-Linux-x86_64.sh


from source code (easy way)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    git clone https://github.com/Amber-MD/pytraj
    cd pytraj
    python ./setup.py install

    # note: pytraj will automatically install cpptraj first.

Depend on your machine, the fresh installation (``libcpptraj`` + ``pytraj``) could take 2 to 4 minutes.

from source code (hard way - expert only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you want to install `libcpptraj.so` by yourself.

- First, download cpptraj::

    git clone https://github.com/Amber-MD/cpptraj
    cd cpptraj
    bash configure -shared -openmp gnu
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

from `pip`
~~~~~~~~~

Since pytraj depends on libcpptraj, install via pip is not easier than intall from source code.

However, if you still want to do it, make sure to install libcpptraj by yourself and set CPPPTRAJHOME (see above step), then::

    pip install https://github.com/Amber-MD/pytraj/archive/master.zip


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

undefined symbol error
~~~~~~~~~~~~~~~~~~~~~~

For end users, install pytraj is very straigh forward by ``python setup.y install``. But
for whom wants to catch up the development of pytraj, you might get ``undefined symbol
error`` when install new code. This happens you need to keep `pytraj` and `cpptraj` syncs. Sometimes `cpptraj` API is changed and you need to
update `pytraj` code (by ``git pull``) and then recompile pytraj from fresh.

If you already tried hard to install but not successful, it's better to remove old pytraj installation (NOT pytraj source
code) and remove all the old `libcpptraj.so` files (come from conda install or from using ``python setup.py install``...)

- First, remove all `libcpptraj.so` files. You can find their dir by using::
      
    locate libcpptraj.so

- Then, remove build directory::

   rm -rf build

- Remove installed folder, example::

   rm -rf /home/anaconda3/lib/python3.4/site-packages/pytraj/

- Recompile `libcpptraj.so`::

  cd cpptraj

  make libcpptraj

- Go back to `pytraj` source::

  python setup.py install

If above steps do not solve your probrem, please contact us.


install ipython and its notebook for interactive data exploration
-----------------------------------------------------------------

`ipython <http://ipython.org/>`_ and its notebook is great program for interactive exloration of MD data.
Curious about how the notebook looks like? check out our `pairwise rmsd tutorial <http://amber-md.github.io/pytraj/doc/build/html/tutorials/tut_pairwise_rmsd.html>`_

If you are using ``anaconda``, just type ``ipython notebook``. If you have not haved ipython and its notebook, try ``conda install ipython``
For further instruction and information about ``ipython-notebook``, please check its website.

install anaconda for all python packages
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
