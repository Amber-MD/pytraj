# adapted slightly from ``mdtraj`` package
# coveralls

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi

if [[ $CPPTRAJ_ANACONDA == "NO" ]]; then
    anaconda -t $TRAVIS_TO_ANACONDA upload --force -u ambermd -p pytraj-dev $HOME/miniconda/conda-bld/linux-64/pytraj-dev-*
fi

# only need to update one version for libcpptraj
# only need to update infrequently
#if [[ "$PYTHON_VERSION" = "3.4" ]]; then
#    anaconda -t $TRAVIS_TO_ANACONDA upload --force -u ambermd -p libcpptraj-dev $HOME/miniconda/conda-bld/linux-64/libcpptraj-dev-*
#fi
