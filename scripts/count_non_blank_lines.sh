#!/bin/sh

cat ./pytraj/*.p* ./examples/*.py ./pytraj/*/*.p* | sed '/^\s*$/d' | wc -l 
