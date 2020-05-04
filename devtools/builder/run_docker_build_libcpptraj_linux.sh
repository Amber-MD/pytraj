#!/usr/bin/env bash
set -e

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=hainm/pytraj-build-box:2019-03

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        ${DOCKER_IMAGE}\
                        bash || exit $?

set -x
set -e
cd /feedstock_root/

if [ ! -d cpptraj ]; then
    python scripts/install_libcpptraj.py github -openmp
else
    # (cd cpptraj && git pull && git clean -fdx .)
    python scripts/install_libcpptraj.py -openmp
fi
EOF
