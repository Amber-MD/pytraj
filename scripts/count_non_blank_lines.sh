#!/bin/sh

cat ./pycpptraj/*.p* ./examples/*.py ./pycpptraj/*/*.p* | sed '/^\s*$/d' | wc -l 
