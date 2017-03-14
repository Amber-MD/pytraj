#!/bin/sh

# Must be in pytraj root folder

# add cpptraj folder
if [ ! -d cpptraj ]; then
    git clone https://github.com/amber-md/cpptraj
else
    (cd cpptraj && git pull && git clean -fdx .)
fi

# Linux build via docker
# sh devtools/ci/run_docker_build_wheels.sh

# Osx build
(cd cpptraj && git clean -fdx .)
export CPPTRAJHOME=`pwd`/cpptraj
python scripts/install_libcpptraj.py
(cd dist && python ../scripts/build_wheel.py pytraj*.tar.gz)
