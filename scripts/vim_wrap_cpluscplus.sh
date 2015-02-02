#!/bin/sh

fname=$1

if [ ! -f $fname ]; then
cat >$fname <<EOF
# distutils: language = c++
EOF
fi

vim $fname
