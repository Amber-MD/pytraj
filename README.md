1. sphinxext was copied from pandas.

2. Install some programs by running:

    $ conda install sphinx numpy matplotlib seaborn numpy jupyter notebook runipy --yes
    $ pip install sphinx_bootstrap_theme

3. Notes

There are some ipython-notebooks that are not automatically run, (such as
interface_with_mdtraj.ipynb) to avoid dependency when building this doc.
Only html file is included.

To run specific notebook without opening it, use `runipy your_notebook.ipynb`
