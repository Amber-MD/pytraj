# distutils: language = c++

from cpp_vector cimport vector as cppvector
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from pytraj.externals.six import string_types

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

    def __iter__(self):
        cdef cppvector[string].const_iterator it
        it = self.thisptr.begin()
        while it != self.thisptr.end():
            yield deref(it)
            incr(it)

    @classmethod
    def copy(cls, ArgList rhs):
        cdef ArgList arglist = ArgList()
        #del arglist.thisptr
        arglist.thisptr = new _ArgList(rhs.thisptr[0])
        return arglist

    def __getitem__(ArgList self, idx):
        return self.thisptr[0][idx]

    def clear_list(self):
        self.thisptr.ClearList()

    def set_list(self, string inputString, char* separator):
        return self.thisptr.SetList(inputString, separator)
    
    def remaining_args(self):
        cdef ArgList remain = ArgList()
        remain.thisptr[0] = self.thisptr.RemainingArgs()
        return remain

    def add_arg(self, string minput):
        self.thisptr.AddArg(minput)

    def list(self):
        return self.thisptr.List()

    @property
    def n_args(self):
        return self.thisptr.Nargs()

    def is_empty(self):
        return self.thisptr.empty()

    def arg_line(self):
        return self.thisptr.ArgLine()

    def mark_arg(self, int arg):
        self.thisptr.MarkArg(arg)

    def check_for_more_args(self):
        return self.thisptr.CheckForMoreArgs()

    def print_list(self):
        self.thisptr.PrintList()

    def print_debug(self):
        self.thisptr.PrintDebug()

    def remove_first_arg(self):
        self.thisptr.RemoveFirstArg()

    def command_is(self, char* cm):
        return self.thisptr.CommandIs(cm)
    
    def get_string_next(self):
        key = self.thisptr.GetStringNext()
        return key.decode()

    def get_string_key(self, c):
        key = self.thisptr.GetStringKey(c.encode())
        return key.decode()

    def get_mask_next(self):
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

    def contains(self, key):
        return self.thisptr.Contains(key.encode())
