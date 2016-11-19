#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=ambermd/manylinux-extra

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

set -e -x
cd /feedstock_root/

for pyver in cp35-cp35m cp35-cp35m cp35-cp35m; do
    export pybin=/opt/python/\${pyver}/bin/
    \$pybin/python -m pip install cython
done

# use python=3.5 for workflow
export pyver=cp35-cp35m
export \$PATH=\$pybin:$PATH

\$pybin/python scripts/install_libcpptraj.py github -openmp
export CPPTRAJHOME=`pwd`/cpptraj/
\$pybin/python devtools/mkrelease
ls .
cd dist
\$pybin/python ../scripts/build_wheels.py pytraj*gz --manylinux-docker
cd ..
EOF
