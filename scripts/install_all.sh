#!/bin/sh

conda install sphinx numpy matplotlib seaborn pandas jupyter notebook runipy --yes
conda install scikit-learn --yes

# osx
# pip install sphinx_bootstrap_theme mpi4py

# linux
conda install mpi4py --yes
pip install sphinx_bootstrap_theme
