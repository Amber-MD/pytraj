from libcpp.vector cimport vector

cdef sshort[:] gdt(double[:] reference, double[:] arr, int conformers, int protlen, int score):
#cdef sshort[:] gdt(pyarray[double] reference, pyarray[double] arr, int conformers, int protlen, int score):
    cdef sshort[:] arr0
    cdef sshort* ptr

    ptr = <sshort*> _gdt(&reference[0], &arr[0], conformers, protlen, score)
    arr0 = <sshort[:conformers]> ptr
    return arr0
