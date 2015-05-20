# distutils: language = c++

from pytraj.cpptraj_dict import DataFormatDict, get_key
from pytraj.externals.six import string_types

cdef class DataFile:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFile()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    @classmethod
    def write_help(cls):
        _DataFile.WriteHelp()

    @classmethod
    def read_options(cls):
        _DataFile.ReadOptions()

    @classmethod
    def write_options(cls):
        _DataFile.WriteOptions()

    @classmethod
    def get_format_from_arg(cls, ArgList a):
        return _DataFile.GetFormatFromArg(a.thisptr[0])

    #@classmethod
    def format_string(self, string t):
        cdef DataFormatType fmt = DataFormatDict[t]
        return self.thisptr.FormatString(fmt)

    def set_precision(self, int widthIn, int precisionIn):
        self.thisptr.SetDataFilePrecision(widthIn, precisionIn)

    def read_data(self, filenameIn, arglist, DataSetList datasetlist):
        return self.thisptr.ReadDataIn(filenameIn.encode(), 
               ArgList(arglist).thisptr[0], datasetlist.thisptr[0])

    def setup_datafile(self, string filenameIn, ArgList argIn, int debugIn):
        return self.thisptr.SetupDatafile(filenameIn, argIn.thisptr[0], debugIn)

    def add_dataset(self, DataSet dataIn):
        return self.thisptr.AddSet(dataIn.baseptr0)

    def removeset(self, DataSet dataIn):
        return self.thisptr.RemoveSet(dataIn.baseptr0)

    def process_args(self, arg):
        cdef string s
        cdef ArgList argIn

        if isinstance(arg, string_types):
            s = <string> arg
            return self.thisptr.ProcessArgs(s)
        elif isinstance(arg, ArgList):
            argIn = <ArgList> arg
            return self.thisptr.ProcessArgs(argIn.thisptr[0])
        else:
            raise ValueError()

    def write_data(self):
        self.thisptr.WriteData()

    def data_set_names(self):
        self.thisptr.DataSetNames()

    def setDFLwrite(self, bint fIn):
        self.thisptr.SetDFLwrite(fIn)

    def DFLwrite(self):
        return self.thisptr.DFLwrite()

    property dtype:
        def __get__(self):
            return get_key(self.thisptr.Type(), DataFormatDict)
