# adapted slightly from ``mdtraj`` package
#coveralls

echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "setup" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


anaconda -t $TRAVIS_TO_ANACONDA upload --force -u ambermd -p pytraj-dev $HOME/miniconda/conda-bld/linux-64/pytraj-dev-*
#binstar -t $TRAVIS_TO_ANACONDA upload --force -u ambermd -p pytraj-dev $HOME/miniconda/conda-bld/linux-64/pytraj-dev-*
