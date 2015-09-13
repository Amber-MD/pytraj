# copied from `ParmEd`
MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b

export PATH=$HOME/miniconda/bin:$PATH
conda install --yes conda-build jinja2 binstar pip anaconda-client
conda config --add channels http://conda.binstar.org/ambermd/

if [ -z "$NO_CYTHON" ]; then
    if [ $PYTHON_VERSION = '3.4' ]; then
        python ./scripts/clean_cythonize_files_git.py
    fi
    conda create -y -n myenv python=$PYTHON_VERSION \
        numpy cython h5py
else
    conda create -y -n myenv python=$PYTHON_VERSION numpy h5py
fi

source activate myenv

# install other packages here
pip install git+git://github.com/ParmEd/ParmEd
pip install memory_profiler
pip install psutil
pip install cclib
