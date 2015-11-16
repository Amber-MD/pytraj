1. sphinxext was copied from pandas.

2. Need

- sphinx: I am using v1.2.3
- numpy
- matplotlib
- seaborn: v0.0.6
- jupyter (conda install jupyter)
- sphinx_bootstrap_theme: I am using v0.4.7
- runipy
- sphinxcontrib-lunrsearch (instant search box)

Most above packages can be installed via conda
    example: conda install matplotlib seaborn numpy jupyter

3. Notes

There are some ipython-notebooks that are not automatically run, (such as
interface_with_mdtraj.ipynb) to avoid dependency when building this doc.
Only html file is included.
