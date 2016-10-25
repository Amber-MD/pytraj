# Release

## Source code
```bash
./devtools/mkrelease

# output file: dist/pytraj-{version}.tar.gz
# e.g: dist/pytraj-1.1.0.tar.gz
```

## Wheel file for pypi

1. build libcpptraj and set CPPTRAJHOME

2. build wheel files
```bash
python ./scripts/build_wheel.py dist/pytraj-{version}.tar.gz
# output: wheelhouse/pytraj-{version}-*.whl
# - upload
twine upload wheelhouse/pytraj-{version}-*.whl
```

## conda
#
Please check [amber-recipes](https://github.com/Amber-MD/amber-recipes/tree/master/pytraj) repo

# See also

[Developer guide](http://amber-md.github.io/pytraj/latest/developer_guide.html)
