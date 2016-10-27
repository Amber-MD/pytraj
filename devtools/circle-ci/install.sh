#!/bin/sh

conda update conda-build -y
conda build devtools/conda-recipe/libcpptraj/ -q
conda build devtools/conda-recipe/pytraj --python=$PYTHON_VERSION -q

conda install `conda build --output devtools/conda-recipe/libcpptraj --python=$PYTHON_VERSION` -y
conda install `conda build --output devtools/conda-recipe/pytraj --python=$PYTHON_VERSION` -y
