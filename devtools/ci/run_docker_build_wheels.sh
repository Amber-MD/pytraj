#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=ambermd/manylinux-extra

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

cd /feedstock_root/
export pyver=cp35-cp35m
export pybin=/opt/python/\${pyver}/bin/
export \$PATH=\$pybin:$PATH
./devtools/mkrelease
cd dist
\$pybin/python ../scripts/build_wheels.py pytraj*gz --manylinux-docker
cd ..
EOF
