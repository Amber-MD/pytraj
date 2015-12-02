.. _conda:

Use conda to manage different Python versions
---------------------------------------------

Supposed you want to switch back and forth among different Python versions (python2.7, 3.4, 3.5, ...), or you
want to try different environments (with or without numpy):  ``conda`` is solution

- install conda: http://conda.pydata.org/docs/installation.html ::

  Note: if you have Anaconda installed, conda is already there.
 
- To create a new enviroment with name of ``your_env_name``, python version=3.4, preinstall ``numpy``, ``scipy`` ::

  conda create -n your_env_name python=3.4 numpy scipy

- activate to newly created enviroment ::

  source activate your_env_name

- if you want to install new program in ``your_env_name``, such as ``matplotlib`` ::

  conda install matplotlib

- return to your original env ::

  source deactivate

Read also: http://conda.pydata.org/docs/using/index.html
