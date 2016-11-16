#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
CONDA=/root/miniconda/bin/conda
DOCKER_IMAGE=centos:5

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

yum -y update
yum -y install gcc \
               patch \
               csh \
               flex \
               wget \
               perl \
               bzip2 \
               libgfortran44.x86_64 \
               make \
               m4 \
               which

# Embarking on 1 case(s).
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b
    export PATH=/root/miniconda/bin:\$PATH
    $CONDA update --yes --all
    $CONDA install --yes conda-build anaconda-client
    $CONDA info
    $CONDA install --yes numpy nomkl zlib bzip2 libnetcdf libgfortran openblas gcc
    pip install auditwheel
    cd /feedstock_root/
    ls .
    sh devtools/ci/test_pip_build.sh
    cd -
    $CONDA build /cpptraj_recipe --quiet || exit 1
    $CONDA build /feedstock_root/devtools/conda-recipe/pytraj --quiet || exit 1
    cp /opt/conda/conda-bld/linux-64/pytraj*  /feedstock_root/
    cp /opt/conda/conda-bld/linux-64/libcpptraj* /feedstock_root/
EOF
