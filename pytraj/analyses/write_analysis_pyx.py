import re
from glob import glob

pxdlist = glob("Analysis_*.pxd")
analysislist = []

for pxd in pxdlist:
    try:
        analysis = re.findall("Analysis_(.+?).pxd", pxd)[0]
    except:
        pass
    analysislist.append(analysis)

exlucdedList = []
for excluded_analysis in exlucdedList:
    analysislist.remove(excluded_analysis)

#print analysislist
text = """# distutils: language = c++
from cython.operator cimport dereference as deref


cdef class Analysis_ANALYSIS_NAME (Analysis):
    def __cinit__(self):
        self.baseptr = <_Analysis*> new _Analysis_ANALYSIS_NAME()
        self.thisptr = <_Analysis_ANALYSIS_NAME*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        \"""return a function-pointer object to be used with AnalysisList class
        \"""
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()
""" 

for analysis in analysislist:
    tmp = text.replace("ANALYSIS_NAME", analysis)
    fname = "Analysis_" + analysis + ".pyx"
    with open(fname, 'w') as fh:
        fh.writelines(tmp)
