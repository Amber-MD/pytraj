Making a Release
----------------

- create git tag
- build binary distribution
```bash
git clone https://github.com/hainm/build_pytraj
cd build_pytraj
sh build_conda.sh
sh build_pip.sh
# conda install anaconda-client
anaconda upload /path/to/conda-installable/binary/file --user ambermd
# pip install twine
twine upload pytraj/dist/wheelhouse/pytraj*.whl
```

Further info
------------
[pytraj website](http://amber-md.github.io/pytraj/latest/developer_guide.html)
