#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=ambermd/pytraj_build
MINICONDA_ROOT=/root/miniconda3
CONDA=$MINICONDA_ROOT/bin/conda

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

export PATH=$MINICONDA_ROOT/bin:\$PATH
$CONDA update --yes --all
cd /feedstock_root/
ls .
sh devtools/ci/test_pip_build.sh
cd -
$CONDA build /cpptraj_recipe --quiet || exit 1
$CONDA build /feedstock_root/devtools/conda-recipe/pytraj --quiet || exit 1
cp $MINICONDA_ROOT/conda-bld/linux-64/pytraj*  /feedstock_root/
cp $MINICONDA_ROOT/conda-bld/linux-64/libcpptraj* /feedstock_root/
EOF
