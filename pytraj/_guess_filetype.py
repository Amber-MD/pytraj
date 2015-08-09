from pytraj.externals.six import PY3


def _guess_filetype(myfile):
    """determine filetype of `myfile` based on its contents
    Use with your own risk
    """
    if PY3:
        mode = 'br'
    else:
        # PY2
        mode = 'r'

    with open(myfile, mode) as fh:
        first_line = fh.read()
        if b'HDF' in first_line[:10]:
            return 'HDF5'
