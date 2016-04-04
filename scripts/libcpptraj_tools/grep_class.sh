#!/bin/sh

cname=$1

grep  $cname $AMBERHOME//AmberTools/src/cpptraj/src/$cname.cpp | sed "/\/\//d"
