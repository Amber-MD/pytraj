# distutils: language = c++

cdef class BaseIOtype:
    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass


cdef class DispatchObject:
    def __cinit__(self):
        self.thisptr = new _DispatchObject()

    def __dealloc__(self):
        del self.thisptr

cdef class FunctPtr:
    def __cinit__(self):
        # just dummy class
        pass

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


cdef class CoordinateInfo:
    def __cinit__(self):
        self.thisptr = new _CoordinateInfo()

    def __dealloc__(self):
        del self.thisptr

    def has_box(self):
        return self.thisptr.HasBox()

    def traj_box(self):
        cdef Box box = Box()
        box.thisptr[0] = self.thisptr.TrajBox()
        return box

    def has_vel(self):
        return self.thisptr.HasVel()

    def has_temp(self):
        return self.thisptr.HasTemp()

    def has_time(self):
        return self.thisptr.HasTime()

    def has_force(self):
        return self.thisptr.HasForce()

    def has_replica_dims(self):
        return self.thisptr.HasReplicaDims()

from pytraj.cpptraj_dict import AccessDict
from pytraj import cpptraj_dict

cdef class CpptrajFile:
    """
    Original cpptraj doc:
    Class to abstract handling of basic file routines.
    """
    def __cinit__(self, *args, **kwd):
        """
        >>> cpptraj = CpptrajFile()
        >>> cpptraj2 = CpptrajFile(cpptraj)
        >>> cpptraj23 = CpptrajFile("test.dat", 'r')
        >>> with CpptrajFile("test2.dat", 'r') as cfile:
        >>>    pass
        """

        cdef CpptrajFile cfile
        if not args and kwd:
            self.thisptr = new _CpptrajFile()
        else:
            if isinstance(args[0], CpptrajFile):
                # ignore other options
                cfile = args[0]
                self.thisptr = new _CpptrajFile(cfile.thisptr[0])
            else:
                self.thisptr = new _CpptrajFile()
                self.open(*args, **kwd)

    def __dealloc__(self):
        """ This is virtual method"""
        #del self.thisptr
        pass

    def __enter__(self):
        return self

    def __exit__(self, arg1, arg2, arg3):
        self.close()

    def open_read(self, filename):
        cdef bint sucess
        cdef int result
        filename = filename.encode()
        self.thisptr.OpenRead(filename)

    def _setup_read(self, filename, int debug):
        filename = filename.encode()
        return self.thisptr.SetupRead(filename, debug)

    def open_write_numbered(self, int numIn):
        return self.thisptr.OpenWriteNumbered(numIn)

    def open_write(self, filename):
        filename = filename.encode()
        return self.thisptr.OpenWrite(filename)

    def _setup_write(self, *args):
        cdef filename
        cdef FileType ftype
        cdef int debug

        if len(args) == 3:
            filename, ftype, debug = args
            filename = filename.encode()
            return self.thisptr.SetupWrite(filename, ftype, debug)
        elif len(args) == 2:
            filename, debug = args
            filename = filename.encode()
            return self.thisptr.SetupWrite(filename, debug)
        else:
            raise ValueError()

    def open(self, filename="", status='r'):
        """just like Python file
        >>> cpptraj = CpptrajFile()
        >>> cpptraj.open("test.dat", 'r')
        """
        status = status.lower()
        filename = filename.encode()

        if status == 'r' or status == 'read': 
            self.thisptr.OpenRead(filename)
        elif status == 'a' or status == 'append':
            self.thisptr.OpenAppend(filename)
        elif status == 'w' or status == 'write':
            self.thisptr.OpenWrite(filename)
        else:
            raise ValueError("wrong open status")

    def open_append(self, filename):
        filename = filename.encode()
        return self.thisptr.OpenAppend(filename)

    def _setup_append(self, filename, int debug):
        filename = filename.encode()
        return self.thisptr.SetupAppend(filename, debug)

    def openfile(self, *args):
        # TODO : what this method does?
        cdef AccessType accessIn
        if not args:
            return self.thisptr.OpenFile()
        else:
            accessIn =args[0]
            return self.thisptr.OpenFile(accessIn)
        
    def close(self):
        self.thisptr.CloseFile()

    def get_line(self):
        return self.thisptr.GetLine()

    def nextline(self):
        # return char*
        return self.thisptr.NextLine()

    @property
    def mode(self):
        key = cpptraj_dict.get_key(self.thisptr.Access(), AccessDict)
        return key.lower()

    def compression(self):
        return self.thisptr.Compression()

    def is_open(self):
        return self.thisptr.IsOpen()

    def filename(self):
        cdef FileName filename = FileName()
        filename.thisptr[0] = self.thisptr.Filename()
        return filename

    def file_size(self):
        return self.thisptr.FileSize()

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

    def __str__(self):
        return self.truncated_name

    def __repr__(self):
        txt = " <atom type: %s>" % self
        return txt 

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

    def __getitem__(self, int idx):
        return self.thisptr.opr_idx(idx)

    @property
    def truncated_name(self):
        """return string"""
        return self.thisptr.Truncated().decode()

    def replace_asterisk(self):
        self.thisptr.ReplaceAsterisk()

    def __str__(self):
        return (self.thisptr.opr_star()).decode()
