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

    @property
    def fullname(self):
        return self.thisptr.Full()
