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
# install stable version?
pip install conda
# conda update conda --yes
conda install --yes conda-build jinja2 anaconda-client pip
#conda config --add channels http://conda.binstar.org/ambermd/

conda create -y -n myenv python=$PYTHON_VERSION numpy cython h5py

source activate myenv
conda install --yes anaconda-client coverage python-coveralls nose
conda install --yes scipy mpi4py libnetcdf

# install other packages here
pip install git+git://github.com/ParmEd/ParmEd
pip install memory_profiler
pip install psutil
pip install cclib
