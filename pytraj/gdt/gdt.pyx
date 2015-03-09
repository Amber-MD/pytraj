from libcpp.vector cimport vector

ctypedef signed short sshort
cdef extern from "./src/gdt_.h":
    sshort * _gdt "gdtCPUOneReference"(double *reference, double *arr,  int conformers, int protlen,int score)


def gdt(double[:] reference, double[:] arr, int conformers, int protlen, int score):
    """gdt(double[:] reference, double[:] arr, int conformers, int protlen, int score)"""
    cdef sshort[:] arr0
    cdef sshort* ptr

    ptr = <sshort*> _gdt(&reference[0], &arr[0], conformers, protlen, score)
    arr0 = <sshort[:conformers]> ptr
    return arr0
