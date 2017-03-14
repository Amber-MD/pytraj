#!/bin/sh

# Must be in pytraj root folder

# add cpptraj folder
python -c "import auditwheel" || exit 1

if [ ! -d cpptraj ]; then
    git clone https://github.com/amber-md/cpptraj
else
    (cd cpptraj && git pull && git clean -fdx .)
fi

# PIP BUILD
#     LINUX
sh devtools/ci/run_docker_build_wheels.sh

#     OSX
(cd cpptraj && git clean -fdx .)
export CPPTRAJHOME=`pwd`/cpptraj
python scripts/install_libcpptraj.py
(cd dist && python ../scripts/build_wheel.py pytraj*.tar.gz)

# CONDA BUILD
#     LINUX
sh devtools/ci/run_docker_build_conda.sh

#     OSX
for pyver in 2.7 3.4 3.5; do
    conda build devtools/conda-recipe/pytraj --py $pyver
    tarfile=`conda build devtools/conda-recipe/pytraj --py $pyver --output`

    build_dir=dist/conda/osx-64
    if [ ! -d $build_dir ]; then
        mkdir -p $build_dir
    fi
    cp $tarfile $build_dir
done
