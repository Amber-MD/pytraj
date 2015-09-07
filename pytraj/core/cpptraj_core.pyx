# distutils: language = c++
from .CpptrajState cimport _CpptrajState, CpptrajState
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarrary
from cpython cimport array
from ..externals.six import string_types
from ..utils import is_array, is_int
from .._cyutils import _int_array1d_like_to_memview
from ..externals.six.moves import range


cdef class AtomMask(object):
    def __cinit__(self, *args):
        cdef int begin_atom, end_atom, atom_num
        cdef string maskstring
        cdef AtomMask rhs_atm
        cdef vector[int] v_int
        cdef int i, max_atoms

        if not args:
            self.thisptr = new _AtomMask()
        elif is_array(args[0]) or isinstance(args[0], (list, tuple, range)):
            self.thisptr = new _AtomMask()
            self.add_selected_indices(args[0])
        else:
            if len(args) == 1:
                if isinstance(args[0], int):
                    atom_num = args[0]
                    self.thisptr = new _AtomMask(atom_num)
                elif isinstance(args[0], AtomMask):
                    rhs_atm = args[0]
                    self.thisptr = new _AtomMask(rhs_atm.thisptr[0])
                elif isinstance(args[0], string_types):
                    # string_types
                    maskstring = args[0].encode("UTF-8")
                    self.thisptr = new _AtomMask(maskstring)
            elif len(args) == 2:
                if is_int(args[0]) and is_int(args[1]):
                    begin_atom, end_atom = args
                    self.thisptr = new _AtomMask(begin_atom, end_atom)
                else:
                    # array-like with max_atoms
                    for i in args[0]:
                        v_int.push_back(i)
                    max_atoms = args[1]
                    self.thisptr = new _AtomMask(v_int, max_atoms)
            else:
                # TODO: better Error
                raise NotImplementedError()

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr

    def selected_indices(self):
        return self.indices

    @property
    def indices(self):
        import numpy as np
        return np.asarray(self.thisptr.Selected())

    @property
    def _indices_view(self):
        cdef vector[int] v = self.thisptr.Selected()
        cdef pyarrary a_empty = pyarrary('i', [])
        cdef int size = v.size()
        cdef pyarrary arr0 = array.clone(a_empty, size, zero=True) 
        cdef int[:] myview = arr0
        cdef int i 

        for i in range(size):
            myview[i] = v[i]

        return myview

    def __iter__(self):
        cdef cppvector[int].const_iterator it
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            yield deref(it)
            incr(it)

    @property
    def n_selected(self):
        return self.thisptr.Nselected()

    property n_atoms:
        def __get__(self):
            """the number of selected atoms based on mask"""
            return self.thisptr.Nselected()

    def __getitem__(self, int idx):
        return self.thisptr.index_opr(idx)

    property mask_string:
        def __get__(self):
            return self.thisptr.MaskString()
        def __set__(self, value):
            self.thisptr.SetMaskString(value.encode())

    def mask_expression(self):
        cdef string t
        t = self.thisptr.MaskExpression()
        return t.decode()

    def is_empty(self):
        return self.thisptr.None()

    def invert_mask(self):
        self.thisptr.InvertMask()

    def add_selected_indices(self, arr0):
        """add atom index without sorting

        See Also
        --------
        add_atom
        """
        cdef int[:] int_view
        cdef int i

        # try casting to memview
        if not is_array(arr0):
            int_view = _int_array1d_like_to_memview(arr0)
        else:
            try:
                int_view = arr0
            except:
                # numpy compat
                int_view = arr0.astype('i4')

        for i in range(int_view.shape[0]):
            self.thisptr.AddSelectedAtom(int_view[i])

    def add_atom(self,int atom_num):
        """add atom index and sort"""
        self.thisptr.AddAtom(atom_num)

    def add_atoms(self, vector[int] v):
        self.thisptr.AddAtoms(v)

    def add_atom_range(self, int begin, int end):
        self.thisptr.AddAtomRange(begin, end)


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


cdef extern from "Command.h": 
    ctypedef enum RetType "Command::RetType":
        pass
    cdef cppclass _Command "Command":
        @staticmethod
        RetType ProcessInput(_CpptrajState&, const string&)

cdef class Command:
    cdef _Command* thisptr

    def __cinit__(self):
        self.thisptr = new _Command()

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def get_state(cls, trajin_text):
        cdef CpptrajState cppstate = CpptrajState()
        trajin_text = trajin_text.encode()
        _Command.ProcessInput(cppstate.thisptr[0], trajin_text)
        return cppstate
