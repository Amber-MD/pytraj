# distutils: language = c++
from libcpp.string cimport string
from pytraj.BaseIOtype cimport *
from pytraj.ArgList cimport *


cdef extern from "FileTypes.h": 
    ctypedef int FileFormatType "FileTypes::FileFormatType"

    ctypedef struct AllocToken "FileTypes::AllocToken":
        const char* Description
        _BaseIOtype.HelpType ReadHelp
        _BaseIOtype.HelpType WriteHelp
        _BaseIOtype.AllocatorType Alloc

    ctypedef struct KeyToken "FileTypes::KeyToken":
        FileFormatType Type
        const char* Key
        const char* Extension

    cdef cppclass _FileTypes "FileTypes":
        FileFormatType GetFormatFromArg(KeyPtr, _ArgList&, FileFormatType)
        FileFormatType GetFormatFromString(KeyPtr, const string&, FileFormatType)
        string GetExtensionForType(KeyPtr, FileFormatType)
        FileFormatType GetTypeFromExtension(KeyPtr, const string&, FileFormatType)
        const char * FormatDescription(AllocPtr, FileFormatType)
        _BaseIOtype * AllocIO(AllocPtr, FileFormatType, bint)
        void ReadOptions(KeyPtr, AllocPtr, FileFormatType)
        void WriteOptions(KeyPtr, AllocPtr, FileFormatType)

cdef class FileTypes:
    cdef _FileTypes* thisptr
