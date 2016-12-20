#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=hainm/pytraj-manylinux-build-box

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        ${DOCKER_IMAGE}\
                        bash || exit $?

set -x
cd /feedstock_root/

export python=/opt/python/cp35-cp35m/bin/python
\$python scripts/install_libcpptraj.py github -openmp
\$python devtools/mkrelease
cd dist
\$python ../scripts/build_wheel.py pytraj*gz --manylinux-docker
cd ..
EOF
