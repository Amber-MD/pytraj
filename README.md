1. sphinxext was copied from pandas.

2. Install some programs by running:

    $ conda install sphinx numpy matplotlib seaborn jupyter notebook runipy --yes
    $ pip install sphinx_bootstrap_theme line_profiler memory_profiler

3. Notes

There are some ipython-notebooks that are not automatically run, (such as
interface_with_mdtraj.ipynb) to avoid dependency when building this doc.
Only html file is included.

To run specific notebook without opening it, use `runipy your_notebook.ipynb`

4. If you want to include html file in .rst file, please make a copy to latest/ folder too
(have not figured out why sphinx did not properly copy files yet)
