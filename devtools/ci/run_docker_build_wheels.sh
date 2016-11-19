#!/usr/bin/env bash

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
DOCKER_IMAGE=ambermd/manylinux-extra

docker info
cat << EOF | docker run -i \
                        -v ${FEEDSTOCK_ROOT}:/feedstock_root \
                        -a stdin -a stdout -a stderr \
                        $DOCKER_IMAGE \
                        bash || exit $?

set -x
cd /feedstock_root/

for pyver in cp27-cp27m cp34-cp34m cp35-cp35m; do
    export pybin=/opt/python/\${pyver}/bin/
    \$pybin/python -m pip install pip --upgrade
    \$pybin/python -m pip install cython
    \$pybin/python -m pip install numpy
done

# use python=3.5 for workflow
export pyver=cp35-cp35m
export \$PATH=\$pybin:$PATH

\$pybin/python scripts/install_libcpptraj.py github -openmp
export CPPTRAJHOME=\`pwd\`/cpptraj/
\$pybin/python devtools/mkrelease
cd dist
\$pybin/python ../scripts/build_wheel.py pytraj*gz --manylinux-docker
cd ..
EOF
