# distutils: language = c++
from pytraj.externals.six import string_types


cdef class NameType:
    def __cinit__(self, *args):
        cdef string s
        cdef NameType rhs

        if not args:
            self.thisptr = new _NameType()
        elif len(args) == 1:
            if isinstance(args[0], string_types):
                s = args[0].encode()
                self.thisptr = new _NameType(s)
            elif isinstance(args[0], NameType):
                rhs = args[0]
                self.thisptr = new _NameType(rhs.thisptr[0])
            else:
                raise ValueError()
        else:
            raise ValueError()

    def copy(self):
        cdef NameType nt = NameType()
        del nt.thisptr 
        nt.thisptr = new _NameType(self.thisptr[0])
        return nt

    def __dealloc__(self):
        del self.thisptr

    def to_buffer(self, char* c):
        # TODO : what does this method do?
        self.thisptr.ToBuffer(c)

    def match(self, NameType nt):
        """return bool"""
        return self.thisptr.Match(nt.thisptr[00])

    def __richcmp__(NameType self, arg, int opt):
        # better way?
        cdef char* c
        cdef bytes py_bytes 
        cdef NameType rhs

        if isinstance(arg, string_types):
            py_bytes = arg.encode()
            c = py_bytes
            if opt == 3:
                # != operator
                return self.thisptr[0] != c
            elif opt == 2:
                # == operator
                return self.thisptr[0] == c
        if isinstance(arg, NameType):
            rhs = arg
            if opt == 3:
                # != operator
                return self.thisptr[0] != rhs.thisptr[0]
            elif opt == 2:
                # == operator
                return self.thisptr[0] == rhs.thisptr[0]

    #def  char * operator *(self):

    def __getitem__(self, int idx):
        return self.thisptr.opr_idx(idx)

    def truncated(self):
        """return string"""
        return self.thisptr.Truncated()

    def replace_asterisk(self):
        self.thisptr.ReplaceAsterisk()

    def __str__(self):
        return str(self.thisptr.opr_star())
