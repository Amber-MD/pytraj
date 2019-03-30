#!/bin/sh

# Must be in pytraj root folder

# add cpptraj folder
python -c "import auditwheel" || exit 1

function main(){
    # this function will be run in the end of this script
    devtools/mkrelease
    clone_or_update_cpptraj
    pip_linux
    pip_osx
    conda_linux
    conda_osx
}


function clone_or_update_cpptraj(){
    if [ ! -d cpptraj ]; then
        git clone https://github.com/amber-md/cpptraj
    else
        (cd cpptraj && git pull && git clean -fdx .)
    fi
}


function pip_linux(){
    sh devtools/builder/run_docker_build_wheels_linux.sh
}


function pip_osx(){
    (cd cpptraj && git clean -fdx .)
    export CPPTRAJHOME=`pwd`/cpptraj
    python scripts/install_libcpptraj.py
    (cd dist && python ../scripts/build_wheel.py pytraj*.tar.gz)
}


function conda_linux(){
    sh devtools/builder/run_docker_build_conda_linux.sh
}


function conda_osx(){
    for pyver in 2.7 3.7 3.5 3.6; do
        conda build devtools/conda-recipe/pytraj --py $pyver
        tarfile=`conda build devtools/conda-recipe/pytraj --py $pyver --output`
    
        build_dir=dist/conda/osx-64
        if [ ! -d $build_dir ]; then
            mkdir -p $build_dir
        fi
        cp $tarfile $build_dir
    done
}

main "$@"
