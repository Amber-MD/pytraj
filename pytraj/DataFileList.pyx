# distutils: language = c++
from cython.operator cimport dereference as deref


cdef class DataFileList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFileList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def remove_datafile(self, DataFile dfIn):
        cdef DataFile dfile = DataFile()
        dfile.thisptr[0] = deref(self.thisptr.RemoveDataFile(dfIn.thisptr))
        return dfile

    def remove_data_set(self,DataSet dsIn):
        self.thisptr.RemoveDataSet(dsIn.baseptr0)
        
    def get_datafile(self, datafilename):
        datafilename = datafilename.encode()
        cdef DataFile dfile = DataFile()
        dfile.py_free_mem = False
        dfile.thisptr = self.thisptr.GetDataFile(datafilename)
        return dfile

    def add_datafile(self, datafilename, *args):
        datafilename = datafilename.encode()
        cdef DataFile dfile = DataFile()
        cdef ArgList argIn

        if not args:
            dfile.thisptr[0] = deref(self.thisptr.AddDataFile(datafilename))
        else:
            argIn = args[0]
            dfile.thisptr[0] = deref(self.thisptr.AddDataFile(datafilename, argIn.thisptr[0]))
        return dfile

    def add_dataset(self, datafilename, ds):
        datafilename = datafilename.encode()
        """we need to use alloc method to cast to DataSet object"""
        dfile = self._add_dataset(datafilename, ds.alloc())
        return dfile

    def _add_dataset(self, datafilename, DataSet dsetIn):
        cdef DataFile dfile = DataFile()
        dfile.thisptr[0] = deref(self.thisptr.AddSetToFile(datafilename, dsetIn.baseptr0))
        return dfile

    def listinfo(self):
        self.thisptr.List()

    def write_all_datafiles(self):
        # perhaps pytraj only uses this method
        self.thisptr.WriteAllDF()

    def reset_write_status(self):
        self.thisptr.ResetWriteStatus()

    def process_data_file_args(self, ArgList argIn):
        return self.thisptr.ProcessDataFileArgs(argIn.thisptr[0])
