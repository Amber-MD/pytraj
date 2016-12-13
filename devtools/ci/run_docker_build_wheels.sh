#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=ambermd/amber-build-box

docker info
cat << EOF | docker run -i \
                        --rm \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        ${DOCKER_IMAGE}\
                        bash || exit $?

    set -x
    cd /feedstock_root/

    # using Python from Miniconda
    pip install auditwheel
    which python
    python scripts/install_libcpptraj.py github -openmp
    export CPPTRAJHOME=\`pwd\`/cpptraj/
    python devtools/mkrelease
    cd dist
    python ../scripts/build_wheel.py pytraj*gz
    cd ..
EOF
