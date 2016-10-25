# Release

## Source code
```bash
./devtools/mkrelease

# output file: dist/pytraj-{version}.tar.gz
# e.g: dist/pytraj-1.1.0.tar.gz
```

## Wheel file for pypi
#
```
# - build libcpptraj and set CPPTRAJHOME

# - build wheel files
python ./scripts/build_wheel.py dist/pytraj-{version}.tar.gz

# output: wheelhouse/pytraj-{version}-*.whl

# - upload
twine upload wheelhouse/pytraj-{version}-*.whl
```

## conda
Please check [amber-receipts](https://github.com/Amber-MD/amber-recipes/tree/master/pytraj) repo
