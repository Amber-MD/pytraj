from libcpp.vector cimport vector
from pytraj.utils.check_and_assert import _import_numpy

has_numpy, np = _import_numpy()

ctypedef signed short sshort
cdef extern from "./src/gdt_.h":
    sshort * _gdt "gdtCPUOneReference"(float *reference, float *arr,  int conformers,int protlen,int score)


def gdt(float[:] reference, float[:] arr, int conformers, int protlen, int score):
    """gdt(float[:] reference, float[:] arr, int conformers, int protlen, int score)"""
    cdef sshort[:] arr0
    cdef sshort* ptr

    ptr = <sshort*> _gdt(&reference[0], &arr[0], conformers, protlen, score)
    arr0 = <sshort[:conformers]> ptr
    if has_numpy:
        return np.asarray(arr0)
    else:
        return arr0
