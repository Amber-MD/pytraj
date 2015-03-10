from cpython.array cimport array as pyarray

ctypedef signed short sshort
cdef extern from "./src/gdt_.h":
    sshort * _gdt "gdtCPUOneReference"(double *reference, double *arr,  int conformers, int protlen,int score)

#cdef sshort[:] gdt(pyarray[double] reference, pyarray[double] arr, int conformers, int protlen, int score)
cdef sshort[:] gdt(double[:] reference, double[:] arr, int conformers, int protlen, int score)
