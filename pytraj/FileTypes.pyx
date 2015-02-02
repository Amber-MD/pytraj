# distutils: language = c++


cdef class FileTypes:
    def __cinit__(self):
        self.thisptr = new _FileTypes()

    def __dealloc__(self):
        del self.thisptr

    #def FileFormatType GetFormatFromArg(self,KeyPtr, ArgList, FileFormatType):

    #def FileFormatType GetFormatFromString(self,KeyPtr, string, FileFormatType):

    #def string GetExtensionForType(self,KeyPtr, FileFormatType):

    #def FileFormatType GetTypeFromExtension(self,KeyPtr, string, FileFormatType):

    #def char * FormatDescription(self,AllocPtr, FileFormatType):

    #def BaseIOtype * AllocIO(self,AllocPtr, FileFormatType, bint):

    #def void ReadOptions(self,KeyPtr, AllocPtr, FileFormatType):

    #def void WriteOptions(self,KeyPtr, AllocPtr, FileFormatType):


