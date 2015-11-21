from __future__ import print_function
import sys

'''
python write_analysis_pxd.py {'Analysis_PhiPsi', 'Analysis_TI', 'Analysis_Wavelet'}

'''

template = '''
cdef class {0}(Analysis):
    def __cinit__(self):
        self.baseptr = <_Analysis*> new _{0}()
        self.thisptr = <_{0}*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with AnalysisList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func

    def help(self):
        self.thisptr.Help()
'''

words = ' '.join(sys.argv[1:]).strip("'").strip('{').strip('}').replace(',', ' ').split()

for w in words:
    print(template.format(w))
