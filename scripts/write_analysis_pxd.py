from __future__ import print_function
import sys

'''
python write_analysis_pxd.py {'Analysis_PhiPsi', 'Analysis_TI', 'Analysis_Wavelet'}

'''

template = '''
cdef extern from "{0}.h":
    cdef cppclass _{0} "{0}" (_Analysis) nogil:
        _{0}()
        _DispatchObject * Alloc()
        void Help()


cdef class {0} (Analysis):
    cdef _{0}* thisptr
'''

words = ' '.join(sys.argv[1:]).strip("'").strip('{').strip('}').replace(',', ' ').split()

for w in words:
    print(template.format(w))
