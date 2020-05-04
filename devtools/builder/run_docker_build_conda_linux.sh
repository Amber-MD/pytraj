#!/usr/bin/env bash
set -e

FEEDSTOCK_ROOT=$(cd "$(dirname "$0")/../../"; pwd;)
# DOCKER_IMAGE=hainm/pytraj-build-box:2020-04-24  # centos-7
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

export PATH=\$HOME/miniconda3/bin:\$PATH

for pyver in  3.5 3.6 3.7 3.8; do
    conda build devtools/conda-recipe/pytraj --py \$pyver
    tarfile=\`conda build devtools/conda-recipe/pytraj --py \$pyver --output\`
    echo "\$tarfile"

    build_dir=dist/conda/linux-64
    if [ ! -d \$build_dir ]; then
        mkdir -p \$build_dir
    fi
    cp \$tarfile \$build_dir
done
EOF
