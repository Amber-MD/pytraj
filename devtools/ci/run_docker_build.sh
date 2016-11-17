#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)

# DOCKER_IMAGE=ambermd/pytraj_build
# MINICONDA_ROOT=/root/miniconda3

DOCKER_IMAGE=condaforge/linux-anvil
MINICONDA_ROOT=/opt/conda/

CONDA=$MINICONDA_ROOT/bin/conda
CPPTRAJ_RECIPE=$FEEDSTOCK_ROOT/devtools/conda-recipe/libcpptraj

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -v ${CPPTRAJ_RECIPE}:/cpptraj_recipe \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

export PATH=$MINICONDA_ROOT/bin:\$PATH
export LD_LIBRARY_PATH=$MINICONDA_ROOT/lib:\$LD_LIBRARY_PATH
$CONDA update --yes --all
$CONDA install -y  \
    zlib \
    bzip2 \
    libnetcdf \
    openblas \
    libgfortran \
    cython \
    gcc
$MINICONDA_ROOT/bin/pip install auditwheel
cd /feedstock_root/
ls .
# turn off pip test for now. Got segmentation fault.
sh devtools/ci/test_pip_build.sh
cd -
$CONDA build /cpptraj_recipe --quiet || exit 1
$CONDA build /feedstock_root/devtools/conda-recipe/pytraj_recipe --quiet || exit 1
cp $MINICONDA_ROOT/conda-bld/linux-64/pytraj*  /feedstock_root/
cp $MINICONDA_ROOT/conda-bld/linux-64/libcpptraj* /feedstock_root/
EOF
