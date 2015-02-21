#!/bin/sh

cat ./pytraj/*.p* ./examples/*.py ./pytraj/*/*.p* ./tests/*.py | sed '/^\s*$/d' | wc -l 
