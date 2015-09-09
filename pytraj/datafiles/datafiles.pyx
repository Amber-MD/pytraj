# distutils: language = c++

from pytraj.cpptraj_dict import DataFormatDict, get_key
from pytraj.externals.six import string_types
from cython.operator cimport dereference as deref


cdef class DataFile:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFile()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def format_string(self, string t):
        cdef DataFormatType format = DataFormatDict[t]
        return self.thisptr.FormatString(format)

    def set_precision(self, int widthIn, int precisionIn):
        self.thisptr.SetDataFilePrecision(widthIn, precisionIn)

    def read_data(self, filenameIn, arglist, DatasetList datasetlist):
        return self.thisptr.ReadDataIn(filenameIn.encode(), 
               ArgList(arglist).thisptr[0], datasetlist.thisptr[0])

    def add_dataset(self, Dataset dataIn):
        return self.thisptr.AddSet(dataIn.baseptr0)

    def write_data(self):
        self.thisptr.WriteData()

    property dtype:
        def __get__(self):
            return get_key(self.thisptr.Type(), DataFormatDict)

cdef class DataFileList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFileList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def get_datafile(self, datafilename):
        datafilename = datafilename.encode()
        cdef DataFile dfile = DataFile()
        dfile.py_free_mem = False
        dfile.thisptr = self.thisptr.GetDataFile(datafilename)
        return dfile

    def add_datafile(self, datafilename, *args):
        """
        Parameters
        ----------
        datafilename : str
            output's name
        args : ArgList object, optional
        """
        datafilename = datafilename.encode()
        cdef DataFile dfile = DataFile()
        cdef ArgList argIn

        if not args:
            dfile.thisptr[0] = deref(self.thisptr.AddDataFile(datafilename))
        else:
            argIn = args[0]
            dfile.thisptr[0] = deref(self.thisptr.AddDataFile(datafilename, argIn.thisptr[0]))
        dfile.py_free_mem = False # let DataFileList free memory
        return dfile

    def add_dataset(self, datafilename, Dataset dsetIn):
        cdef DataFile dfile = DataFile()
        datafilename = datafilename.encode()
        dfile.thisptr = self.thisptr.AddSetToFile(datafilename, dsetIn.baseptr0)
        dfile.py_free_mem = False
        return dfile

    def write_all_datafiles(self):
        # perhaps pytraj only uses this method
        self.thisptr.WriteAllDF()
