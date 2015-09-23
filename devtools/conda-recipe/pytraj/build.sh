#!/bin/sh

export CPPTRAJHOME=$PREFIX/lib/cpptraj/
cp -r $RECIPE_DIR/../../.. $SRC_DIR
$PYTHON setup.py clean
$PYTHON setup.py install openmp
