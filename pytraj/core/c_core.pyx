# distutils: language = c++
import numpy as np
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cpython.array cimport array as pyarray
from cpython cimport array

from ..utils import is_array, is_int
from ..cyutils import _int_array1d_like_to_memview
from ..externals.six import string_types
from ..externals.six.moves import range
from pytraj.externals.six import string_types
from pytraj import c_dict

__all__ = ['AtomMask',
           'FileName',
           'CpptrajFile', 'NameType', 'Command',
           'CpptrajState', 'ArgList', ]

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

    def add_atom(self, int atom_num):
        """add atom index and sort"""
        self.thisptr.AddAtom(atom_num)

    def add_atoms(self, vector[int] v):
        self.thisptr.AddAtoms(v)

    def add_atom_range(self, int begin, int end):
        self.thisptr.AddAtomRange(begin, end)


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


cdef class Command:
    cdef _Command* thisptr

    def __cinit__(self):
        self.thisptr = new _Command()

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def __enter__(self):
        _Command.Init()
        return self

    @classmethod
    def __exit__(self,*args):
        _Command.Free()

    @classmethod
    def get_state(cls, trajin_text):
        cdef CpptrajState cppstate = CpptrajState()
        trajin_text = trajin_text.encode()
        _Command.ProcessInput(cppstate.thisptr[0], trajin_text)
        return cppstate

    @classmethod
    def get_state_from_string(cls, txt):
        cdef CpptrajState cppstate = CpptrajState()
        trajin_text = txt.encode()
        _Command.Dispatch(cppstate.thisptr[0], trajin_text)
        return cppstate

    @classmethod
    def dispatch(cls, CpptrajState state, line):
        _Command.Dispatch(state.thisptr[0], line.encode())


cdef class CpptrajState:
    """
    CpptrajState hold all data per cpptraj run. This class is for internal use.
    Check example

    Examples
    --------
    >>> import pytraj as pt
    >>> text = '''
    parm tz2.parm7
    trajin tz2.nc
    rms
    distance :2 :3
    '''
    >>> state = pt.load_cpptraj_state(text)
    >>> state.run()
    CpptrajState, include:
    <datasetlist: 3 datasets>
    >>> print(state.data[1])
    <pytraj.datasets.DatasetDouble: size=101, key=RMSD_00001>
    values:
    [  2.43182129e-07   4.01623189e+00   6.41421043e+00 ...,   8.27504991e+00
       8.19405473e+00   7.77917637e+00]
    >>> print(state.data.keys())
    ['tz2.parm7', 'RMSD_00001', 'Dis_00002']

    """

    def __cinit__(self):
        self.thisptr = new _CpptrajState()
        self.datafilelist = DataFileList(_own_memory=False)
        self.datasetlist = DatasetList(_own_memory=False)

        # cpptraj will take care of memory deallocating from self.thisptr.PFL(FL, DSL, DFL)
        # We don't free memory again
        # (example: self.toplist.thisptr and self.thisptr.PFL() point to the same address)
        # create memory view
        self.datasetlist.thisptr = &self.thisptr.DSL()
        self.datafilelist.thisptr = &self.thisptr.DFL()

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr

    def __str__(self):
        return 'CpptrajState, include:\n' + \
            '<datasetlist: {0} datasets>'.format(len(self.data))

    property data:
        def __get__(self):
            return self.datasetlist

        def __set__(self, DatasetList dslist):
            # do not del memory. don't know why got double-free mem if doing so.
            #del self.datasetlist.thisptr
            self.datasetlist.thisptr = dslist.thisptr

    def __repr__(self):
        return str(self)

    def compute(self):
        self.thisptr.Run()
        return self

    def run(self):
        self.thisptr.Run()
        return self

def _load_batch(txt, traj=None):
    '''return CpptrajState.

    txt can be text or list of command strings
    '''
    cdef CpptrajState state
    state = CpptrajState()

    if isinstance(txt, (list, tuple)):
        lines = txt
    else:
        lines = [line.lstrip().rstrip() for line in txt.split('\n') if line.strip() != '']

    if traj is not None:
        lines_0 = ['parm %s' % traj.top.filename]

        for fname, frame_slice in zip(traj.filelist, traj._frame_slice_list):
            if len(frame_slice) == 3:
                start, stop, step = frame_slice
            elif len(frame_slice) == 2:
                start, stop = frame_slice
                step = 1
            else:
                raise ValueError('invalid frame_slice')
            if stop == -1:
                _stop = 'last'
            elif stop < -1:
                raise RuntimeError(
                    'does not support negative stop for load_batch (except -1 (last))')
            else:
                _stop = stop
            # add 1 to start since cpptraj ise 1-based index for trajin
            start = start + 1
            lines_0.append(
                'trajin {0} {1} {2} {3}\n'.format(
                    fname,
                    str(start),
                    str(_stop),
                    str(step)))

        # add parm, trajin to lines
        lines = lines_0 + lines
    else:
        lines = lines

    _Command.Init()
    for idx, line in enumerate(lines):
        if not line.startswith('#'):
            _Command.Dispatch(state.thisptr[0], line.encode())
    _Command.Free()
    return state


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
