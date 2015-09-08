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
        cdef DataFormatType format = DataFormatDict[t]
        return self.thisptr.FormatString(format)

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

    def remove_dataset(self,DataSet dsIn):
        self.thisptr.RemoveDataSet(dsIn.baseptr0)
        
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

    def add_dataset(self, datafilename, DataSet dsetIn):
        cdef DataFile dfile = DataFile()
        datafilename = datafilename.encode()
        dfile.thisptr = self.thisptr.AddSetToFile(datafilename, dsetIn.baseptr0)
        dfile.py_free_mem = False
        return dfile

    def info(self):
        from pytraj import set_world_silent
        set_world_silent(False)
        self.thisptr.List()
        set_world_silent(True)

    def write_all_datafiles(self):
        # perhaps pytraj only uses this method
        self.thisptr.WriteAllDF()

    def reset_write_status(self):
        self.thisptr.ResetWriteStatus()

    def process_data_file_args(self, ArgList argIn):
        return self.thisptr.ProcessDataFileArgs(argIn.thisptr[0])
