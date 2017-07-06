import os
import sys
import shutil
from contextlib import contextmanager


@contextmanager
def temporarily_move_libcpptraj(libcpptraj):
    old_dir = os.path.dirname(libcpptraj)
    new_dir = 'tmp'
    ext = '.so' if not sys.platform.startswith('darwin') else '.dylib'
    try:
        os.mkdir('./tmp')
    except OSError:
        pass
    shutil.move(libcpptraj, new_dir)
    yield
    shutil.move(new_dir + '/libcpptraj' + ext, old_dir)
    shutil.rmtree(new_dir)


if __name__ == '__main__':
    libcpptraj = '../cpptraj/lib/libcpptraj.dylib'
    with temporarily_move_libcpptraj(libcpptraj):
        print("hello")
