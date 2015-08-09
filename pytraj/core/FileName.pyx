# distutils: language = c++

cdef class FileName:

    def __cinit__(self, myname=''):
        self.thisptr = new _FileName()
        myname = myname.encode()
        self.thisptr.SetFileName(myname)

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        filename = self.fullname
        filename = filename.decode()
        return filename

    def __repr__(self):
        return self.__str__()

    def set_filename(self, mystring, ext=False):
        mystring = mystring.encode()
        if ext:
            self.thisptr.SetFileNameWithExpansion(mystring)
        else:
            self.thisptr.SetFileName(mystring)

    def clear(self):
        self.thisptr.clear()

    def match_full_or_base(self, mystring):
        mystring = mystring.encode()

        return self.thisptr.MatchFullOrBase(mystring)

    @property
    def fullname(self):
        return self.thisptr.Full()

    @property
    def base(self):
        return self.thisptr.Base()

    property ext:
        def __get__(self):
            return self.thisptr.Ext()

    def is_empty(self):
        return self.thisptr.empty()

    @property
    def compress(self):
        return self.thisptr.Compress()

    def dir_prefix(self):
        return self.thisptr.DirPrefix()
