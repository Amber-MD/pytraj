# distutils: language = c++
import numpy as np
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray
from cpython cimport array

from ..utils import is_array, is_int
from .._cyutils import _int_array1d_like_to_memview
from ..externals.six import string_types
from ..externals.six.moves import range
from pytraj.externals.six import string_types
from pytraj.cpptraj_dict import AccessDict
from pytraj import cpptraj_dict

__all__ = ['AtomMask', 'BaseIOtype', 'DispatchObject',
           'FunctPtr', 'FileName', 'CoordinateInfo',
           'CpptrajFile', 'NameType', 'Command',
           'CpptrajState', 'ArgList',]

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
        return np.asarray(self.thisptr.Selected())

    @property
    def _indices_view(self):
        cdef vector[int] v = self.thisptr.Selected()
        cdef pyarray a_empty = pyarray('i', [])
        cdef int size = v.size()
        cdef pyarray arr0 = array.clone(a_empty, size, zero=True) 
        cdef int[:] myview = arr0
        cdef int i 

        for i in range(size):
            myview[i] = v[i]

        return myview

    def __iter__(self):
        cdef vector[int].const_iterator it
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

cdef class CpptrajState:
    """
    CpptrajState hold instances of:
    + TopologyList
    + DataSetList (having output data)
    + DataFileList

    """
    def __cinit__(self):
        self.thisptr = new _CpptrajState()
        self.toplist = TopologyList(py_free_mem=False)
        self.datafilelist = DataFileList(py_free_mem=False)
        self.datasetlist = DataSetList(py_free_mem=False)

        # cpptraj will take care of memory deallocating from self.thisptr.PFL(FL, DSL, DFL)
        # We don't free memory again 
        # (example: self.toplist.thisptr and self.thisptr.PFL() point to the same address)
        # create memory view
        self.toplist.thisptr = self.thisptr.PFL()
        self.datasetlist.thisptr = self.thisptr.DSL()
        self.datafilelist.thisptr = self.thisptr.DFL()

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
    
    def is_empty(self):
        return self.thisptr.EmptyState()

    def add_trajin(self, arg_or_filename, is_ensemble=None):
        # TODO: add trajector instance?
        cdef string filename
        cdef ArgList argIn
        
        if is_ensemble is not None:
            # reading ensemble
            if isinstance(arg_or_filename, ArgList):
                argIn = arg_or_filename
            elif isinstance(arg_or_filename, string_types):
                argIn = ArgList(arg_or_filename)
            else:
                raise ValueError("")
            self.thisptr.AddTrajin(argIn.thisptr[0], is_ensemble)
        elif isinstance(arg_or_filename, string_types):
            # reading single file
            filename = arg_or_filename.encode()
            self.thisptr.AddTrajin(filename)
        else:
            raise NotImplementedError()

    def run_analyses(self):
        return self.thisptr.RunAnalyses()

    def add_trajout(self, arg):
        """add trajout file
        
        Parameters
        ---------
        arg : str or ArgList object
        """
        cdef string filename
        cdef ArgList arglist

        if isinstance(arg, ArgList):
            arglist = arg
            return self.thisptr.AddTrajout(arglist.thisptr[0])
        elif isinstance(arg, string_types):
            filename = arg.encode()
            return self.thisptr.AddTrajout(filename)
        else:
            raise NotImplementedError()

    def add_reference(self, *args):
        """
        Parameters
        ---------
        filename : str
        arg : ArgList object, optional
        """
        cdef string name
        cdef ArgList arglist

        if len(args) == 1:
            if isinstance(args[0], string_types):
                name =  args[0].encode(0)
                self.thisptr.AddReference(name)
            else:
                raise NotImplementedError()
        elif len(args) == 2:
                name =  args[0].encode(0)
                if isinstance(args[1], string_types):
                    arglist = ArgList(args[1])
                else:
                    arglist = <ArgList> args[1]
                self.thisptr.AddReference(name, arglist.thisptr[0])
        else:
            raise NotImplementedError()

    def add_action(self, actobj, arglist):
        """
        Parameters
        ---------
        actobj : Action object or str
        arglist : ArgList object or str
        """
        # need to explicit casting to FunctPtr because self.thisptr.AddAction need to know type 
        # of variables
        cdef FunctPtr alloc_funct
        cdef ArgList _arglist 

        if isinstance(actobj, string_types):
            # if actobj is string, make Action object
            # then cast to FunctPtr
            from pytraj.action_dict import ADICT
            alloc_funct = ADICT[actobj]().alloc()
        else:
            alloc_funct = <FunctPtr> actobj.alloc()

        if isinstance(arglist, string_types):
            _arglist = ArgList(arglist)
        elif isinstance(arglist, ArgList):
            _arglist = arglist
        else:
            raise ValueError("must be string or ArgList object")

        return self.thisptr.AddAction(alloc_funct.ptr, _arglist.thisptr[0])

    def add_analysis(self, obj, ArgList arglist):
        """temp doc: add_analysis(self, obj, ArgList arglist)
        obj :: Action or Analysis instance
        """
        cdef ArgList _arglist 
        cdef FunctPtr alloc_funct = <FunctPtr> obj.alloc()

        if isinstance(arglist, string_types):
            _arglist = ArgList(arglist)
        elif isinstance(arglist, ArgList):
            _arglist = arglist
        else:
            raise ValueError("must be string or ArgList object")

        return self.thisptr.AddAnalysis(alloc_funct.ptr, _arglist.thisptr[0])

    def list_all(self, ArgList arglist):
        return self.thisptr.ListAll(arglist.thisptr[0])

    def clear_list(self, arglist='all'):
        return self.thisptr.ClearList(ArgList(arglist).thisptr[0])

    def remove_dataset(self, ArgList alist):
        return self.thisptr.RemoveDataSet(alist.thisptr[0])

    def run(self):
        return self.thisptr.Run()

    def write_all_datafiles(self):
        self.thisptr.MasterDataFileWrite()
# distutils: language = c++

cdef class ArgList:
    """
    Original doc from cpptraj:
    ========================
    Class: ArgList
        Hold a list of string arguments and keeps track of their usage.
        Can be set from an input line using SetList(), with arguments separated 
        by a specified delimiter, or arguments can be added one-by-one with AddArg.
        Arguments can be accessed with the various getX routines,
        where X is specific to certain types, e.g. getNextDouble returns
        the next double, getNextMask returns an atom mask expression (i.e.
        it has :, @, % characters etc). All of the getX routines (along with
        the hasKey routine) mark the argument they access as used, so that
        subsequent calls with these functions will not return the same
        argument over and over. 

    pytraj doc:
    =============
    change cpptraj method's name to python style's name
    (hasKey --> has_key)
    """
    
    def __cinit__(self, *args):
        # TODO: need to read cpptraj code for construtor
        cdef string  sinput
        cdef char* sep
        cdef ArgList rhs

        if not args:
            self.thisptr = new _ArgList()
        else:
            if len(args) == 1:
                if isinstance(args[0], ArgList):
                    rhs = args[0]
                    self.thisptr = new _ArgList(rhs.thisptr[0])
                elif isinstance(args[0], string_types):
                    sinput = args[0].encode("UTF-8")
                    self.thisptr = new _ArgList(sinput)
                else:
                    raise ValueError()
            elif len(args) == 2:
                sinput, sep = args
                self.thisptr = new _ArgList(sinput, sep)
            else:
                raise ValueError()

    def __dealloc__(self):
        del self.thisptr

    @property
    def n_args(self):
        return self.thisptr.Nargs()

    def is_empty(self):
        return self.thisptr.empty()

    def command_is(self, char* cm):
        return self.thisptr.CommandIs(cm)
    
    def get_next_string(self):
        key = self.thisptr.GetStringNext()
        return key.decode()

    def get_string_key(self, c):
        key = self.thisptr.GetStringKey(c.encode())
        return key.decode()

    def get_next_mask(self):
        mask = self.thisptr.GetMaskNext()
        return mask.decode()

    def get_next_tag(self):
        return self.thisptr.getNextTag()

    def get_next_integer(self, defint):
        return self.thisptr.getNextInteger(defint)

    def get_next_double(self, double defd):
        return self.thisptr.getNextDouble(defd)

    def get_key_int(self, char* key, int defint):
        return self.thisptr.getKeyInt(key, defint)

    def get_key_double(self, char* key, double defd):
        return self.thisptr.getKeyDouble(key, defd)

    def has_key(self, char* key):
        return self.thisptr.hasKey(key)
