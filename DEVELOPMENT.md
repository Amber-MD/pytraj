Making a Release
----------------

- create git tag
- build binary distribution
```bash
git clone https://github.com/hainm/build-pytraj
cd build-pytraj
sh build_conda.sh
sh build_pip.sh
# conda install anaconda-client
anaconda upload /path/to/conda-installable/binary/file --user ambermd
# pip install twine
twine upload pytraj/dist/wheelhouse/pytraj*.whl
```

Merge pytraj repo to amber repo
--------------------------------

```bash
# add remote branch, only need to do once
cd $AMBERHOME
git remote add pytraj_github https://github.com/Amber-MD/pytraj.git

# create a new branch so we can merge to master by squashing
git branch merge_pytraj
git checkout merge_pytraj

# merge from github
git fetch pytraj_github
# merge all commmits to single one
git pull -s recursive -X subtree=AmberTools/src/pytraj -X theirs --squash pytraj_github master

# if getting merge conflict
git mergetool
# then accept all changes by "remote" (NOT local)
git commit -m 'msg'

git checkout master
git merge merge_pytraj --squash --no-commit
git commit -m 'put your commit message here'

# if everything went well, push to origin
git push origin master
```

Further info
------------
[pytraj website](http://amber-md.github.io/pytraj/latest/developer_guide.html)
