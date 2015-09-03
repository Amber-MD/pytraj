.. _test_cpptraj_api:

Testing cpptraj API change with pytraj
--------------------------------------

clone ``pytraj`` repo if you have not done yet

.. code-block:: bash

    $ git clone https://github.com/Amber-MD/pytraj
    $ # as usual, make a name for original repo. We often use ``upstream``
    $ cd pytraj
    $ git remote add upstream https://github.com/Amber-MD/pytraj
    $ git remote add drroe_pytraj https://github.com/drroe/pytraj
    $ # create a new branch
    $ git branch for_test_cpptraj_api
    $ git checkout for_test_cpptraj_api
    $ # update travis install to point to your ``cpptraj`` new branch (for example ``new_traj_api``)
    $ # just need to update ./installs/install_cpptraj_git.sh file

Upate ``./installs/install_cpptraj_git.sh`` file ::

    1. change "git clone https://github.com/Amber-MD/cpptraj" to "git clone https://github.com/drroe/cpptraj"
    2. add "git checkout new_traj_api" to that file

It's time to push to github so travis can automatically test (will take about 12 minutes)

.. code-block:: bash
    
    git add ./installs/install_cpptraj_git.sh
    git commit -m 'pytraj: trick to run'
    git push drroe_pytraj for_test_cpptraj_api
