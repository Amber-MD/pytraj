# distutils: language = c++
from __future__ import print_function
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types
from pytraj.utils import is_generator
from pytraj.shared_methods import iterframe_master
from cython.operator cimport dereference as deref

_ALL = globals()

def _get_adict():
    ADICT = {}

    for (cname, cls) in _ALL.items():
        if cname.startswith("Action_"):
            actname = cname.split('Action_')[1]
            # create dict of action class
            ADICT[actname.lower()] = cls
    
    # add some commond words to ADICT
    ADICT['surf_LCPO'] = _ALL['Action_Surf']
    ADICT['surf_lcpo'] = _ALL['Action_Surf']
    ADICT['secstruct'] = _ALL['Action_DSSP']
    ADICT['rms'] = _ALL['Action_Rmsd']
    ADICT['superpose'] = _ALL['Action_Rmsd']
    ADICT['drmsd'] = _ALL['Action_DistRmsd']
    ADICT["lipidorder"] = _ALL['Action_OrderParameter']
    ADICT["rog"] = _ALL['Action_Radgyr']
    ADICT["stfcdiffusion"] = _ALL['Action_STFC_Diffusion']
    ADICT["symmrmsd"] = _ALL['Action_SymmetricRmsd']

    return ADICT

class ActionDict:

    def __init__(self):
        self.adict = _get_adict()
        self.action_holder = None

    def __getitem__(self, key):
        # return Action object
        # why do we need action_holder?
        # should we use dict in command.cpp in cpptraj for mapping keyword
        # ('Action_DSSP' --> secstruct)
        self.action_holder = self.adict[key]()
        return self.action_holder

    def __del__(self):
        del self.action_holder

    def keys(self):
        return sorted(self.adict.keys())

cdef class Action:
    '''interface to Cpptraj's Action. For internal use.

    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')
    >>> from pytraj.c_action.c_action import Action_Radgyr
    >>> from pytraj.datasets import DatasetList as CpptrajDataSetList
    >>> from pytraj.datafiles import DataFileList
    >>> dslist = CpptrajDataSetList()
    >>> act = Action_Radgyr(command='@CA', dslist=dslist, top=traj.top)
    >>> act.compute(traj)
    '''

    def __cinit__(self):
        # don't directly create instance of this ABC class.
        self.n_frames = 0
        self.top_is_processed = False

    def __init__(self, command='', Topology top=None, DatasetList dslist=None,
                 DataFileList dflist=DataFileList()):
        # __init__ will be called after __cinit__
        # create __init__ to avoid segmentation fault (why? not sure why)
        # don't directly create instance of this ABC class.

        self._command = command
        self._dslist = dslist
        self._dflist = dflist

        if top is not None and dslist is not None:
            self.read_input(command, top=top, dslist=dslist, dflist=dflist)
            self.setup(top)

    def __dealloc__(self):
        # should I del pointer here or in subclass?
        #del self.baseptr
        pass

    def __del__(self):
        del self.baseptr

    def __str__(self):
        txt = "<pytraj.actions.CpptrajActions.%s>" % (self.__class__.__name__)
        return txt

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args, **kwd):
        """
        >>> from pytraj import *
        >>> traj = io.load("../tz2.nc", "../tz2.parm7")
        >>> dslist = DatasetList.DatasetList()
        >>> adict['jcoupling']("outfile Jcoupling.dat kfile Karplus.txt", traj[0], traj.top, dslist=dslist)
        """
        return self._master(*args, **kwd)

    @makesureABC("Action")
    def read_input(self, command='',
                   top=Topology(),
                   DatasetList dslist=DatasetList(),
                   DataFileList dflist=DataFileList(),
                   int debug=0):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Action_Box().help())
        top : Topology
        dslist : DatasetList instance, default=DatasetList()
        dflist : DataFileList instance, default=DataFileList()
        debug : int, default=0
            debug option from cpptraj. (Do we need this?)
        """
        cdef ArgList arglist
        cdef RetType i_fail
        cdef _ActionInit actioninit_
        actioninit_ = _ActionInit(dslist.thisptr[0], dflist.thisptr[0])

        self.top = top

        if isinstance(command, string_types):
            #command = command.encode("UTF-8")
            arglist = ArgList(command)
        elif isinstance(command, ArgList):
            arglist = <ArgList> command

        i_fail = self.baseptr.Init(arglist.thisptr[0],
                                   actioninit_,
                                   debug)

        if i_fail != OK:
            # check before compute to avoid segfault
            raise ValueError("")
        else:
            return i_fail

    @makesureABC("Action")
    def setup(self, Topology top=Topology(), crdinfo={}, n_frames_t=0, get_new_top=False):
        """pass coordinate_info

        Parameters:
        ----------
        top : Topology
        coordinate_info : a Python dict
        n_frames_t : number of frames associated with Topology
        """
        self.top_is_processed = True
        cdef _ActionSetup actionsetup_
        cdef CoordinateInfo crdinfo_
        cdef Box box
        cdef bint has_velocity, has_time, has_force
        cdef Topology new_top = Topology()

        crdinfo2 = dict((k, v) for k, v in crdinfo.items())

        if 'box' not in crdinfo2:
            crdinfo2['box'] = top.box

        crdinfo_ = CoordinateInfo(crdinfo2)

        actionsetup_ = _ActionSetup(top.thisptr, crdinfo_.thisptr[0], n_frames_t)
        status = self.baseptr.Setup(actionsetup_)

        if status == ERR:
            # cpptraj have a bunch of options, so we only check if there is
            # ERR
            raise RuntimeError('failed to setup action')

        if get_new_top:
            new_top._own_memory = False
            new_top.thisptr[0] = actionsetup_.Top()
            return new_top

    @makesureABC("Action")
    def compute(self, current_frame=None, get_new_frame=False):
        """Perform action on Frame

        Parameters
        ----------
        current_frame : Frame instance need to be processed, default=Frame()
        itx : int, frame index
        get_new_frame : bool
        """
        # debug
        cdef Frame frame, new_frame
        cdef int i
        cdef object traj
        cdef _ActionFrame actframe_

        if isinstance(current_frame, Frame):
            frame = <Frame> current_frame
            actframe_ = _ActionFrame(frame.thisptr, self.n_frames)
            self.baseptr.DoAction(self.n_frames, actframe_)
            self.n_frames += 1

            if get_new_frame:
                new_frame = Frame()
                new_frame.thisptr[0] = actframe_.ModifyFrm()
                return new_frame
        else:
            for frame in iterframe_master(current_frame):
                self.compute(frame)

    @makesureABC("Action")
    def post_process(self):
        self.baseptr.Print()

    def _master(self, command='',
                current_frame=Frame(),
                top=Topology(),
                dslist=DatasetList(),
                dflist=DataFileList(),
                int debug=0):
        """create shortcut
        """
        self.read_input(command=command,
                        top=top,
                        dslist=dslist,
                        dflist=dflist, debug=debug)

        self.setup(top=top)
        self.compute(current_frame)
        return dslist

    def reset_counter(self):
        self.n_frames = 0


cdef class Action_Angle(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Angle()
        self.thisptr = <_Action_Angle*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Distance(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Distance()
        self.thisptr = <_Action_Distance*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Rmsd(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Rmsd()
        self.thisptr = <_Action_Rmsd*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()

cdef class Action_Align(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Align()
        self.thisptr = <_Action_Align*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Dihedral(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Dihedral()
        self.thisptr = <_Action_Dihedral*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_AtomMap(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomMap()
        self.thisptr = <_Action_AtomMap*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Strip(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Strip()
        self.thisptr = <_Action_Strip*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Unstrip(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Unstrip()
        self.thisptr = <_Action_Unstrip*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_DSSP(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DSSP()
        self.thisptr = <_Action_DSSP*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Center(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Center()
        self.thisptr = <_Action_Center*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Hbond(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Hbond()
        self.thisptr = <_Action_Hbond*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Image(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Image()
        self.thisptr = <_Action_Image*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Surf(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Surf()
        self.thisptr = <_Action_Surf*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Radgyr(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Radgyr()
        self.thisptr = <_Action_Radgyr*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Mask(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Mask()
        self.thisptr = <_Action_Mask*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Closest(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Closest()
        self.thisptr = <_Action_Closest*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_NAstruct(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NAstruct()
        self.thisptr = <_Action_NAstruct*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Pucker(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Pucker()
        self.thisptr = <_Action_Pucker*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Outtraj(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Outtraj()
        self.thisptr = <_Action_Outtraj*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Average(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Average()
        self.thisptr = <_Action_Average*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Radial(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Radial()
        self.thisptr = <_Action_Radial*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_DistRmsd(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DistRmsd()
        self.thisptr = <_Action_DistRmsd*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Jcoupling(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Jcoupling()
        self.thisptr = <_Action_Jcoupling*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Pairwise(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Pairwise()
        self.thisptr = <_Action_Pairwise*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Molsurf(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Molsurf()
        self.thisptr = <_Action_Molsurf*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_CheckStructure(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CheckStructure()
        self.thisptr = <_Action_CheckStructure*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_RunningAvg(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_RunningAvg()
        self.thisptr = <_Action_RunningAvg*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_AtomicFluct(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomicFluct()
        self.thisptr = <_Action_AtomicFluct*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Watershell(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Watershell()
        self.thisptr = <_Action_Watershell*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Contacts(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Contacts()
        self.thisptr = <_Action_Contacts*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Vector(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Vector()
        self.thisptr = <_Action_Vector*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Principal(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Principal()
        self.thisptr = <_Action_Principal*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Matrix(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Matrix()
        self.thisptr = <_Action_Matrix*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_LIE(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_LIE()
        self.thisptr = <_Action_LIE*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Grid(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Grid()
        self.thisptr = <_Action_Grid*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_GridFreeEnergy(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_GridFreeEnergy()
        self.thisptr = <_Action_GridFreeEnergy*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Dipole(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Dipole()
        self.thisptr = <_Action_Dipole*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Projection(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Projection()
        self.thisptr = <_Action_Projection*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_ClusterDihedral(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_ClusterDihedral()
        self.thisptr = <_Action_ClusterDihedral*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Unwrap(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Unwrap()
        self.thisptr = <_Action_Unwrap*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Diffusion(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Diffusion()
        self.thisptr = <_Action_Diffusion*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_DNAionTracker(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DNAionTracker()
        self.thisptr = <_Action_DNAionTracker*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Scale(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Scale()
        self.thisptr = <_Action_Scale*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_RandomizeIons(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_RandomizeIons()
        self.thisptr = <_Action_RandomizeIons*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_AutoImage(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AutoImage()
        self.thisptr = <_Action_AutoImage*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_STFC_Diffusion(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_STFC_Diffusion()
        self.thisptr = <_Action_STFC_Diffusion*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_AtomicCorr(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomicCorr()
        self.thisptr = <_Action_AtomicCorr*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Bounds(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Bounds()
        self.thisptr = <_Action_Bounds*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Rotate(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Rotate()
        self.thisptr = <_Action_Rotate*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Translate(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Translate()
        self.thisptr = <_Action_Translate*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Box(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Box()
        self.thisptr = <_Action_Box*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_CreateCrd(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CreateCrd()
        self.thisptr = <_Action_CreateCrd*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_MultiDihedral(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MultiDihedral()
        self.thisptr = <_Action_MultiDihedral*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_MakeStructure(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MakeStructure()
        self.thisptr = <_Action_MakeStructure*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_SymmetricRmsd(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_SymmetricRmsd()
        self.thisptr = <_Action_SymmetricRmsd*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Volmap(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Volmap()
        self.thisptr = <_Action_Volmap*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Spam(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Spam()
        self.thisptr = <_Action_Spam*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Temperature(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Temperature()
        self.thisptr = <_Action_Temperature*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Gist(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Gist()
        self.thisptr = <_Action_Gist*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Density(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Density()
        self.thisptr = <_Action_Density*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_PairDist(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_PairDist()
        self.thisptr = <_Action_PairDist*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_OrderParameter(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_OrderParameter()
        self.thisptr = <_Action_OrderParameter*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_FixAtomOrder(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_FixAtomOrder()
        self.thisptr = <_Action_FixAtomOrder*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_NMRrst(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NMRrst()
        self.thisptr = <_Action_NMRrst*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_FilterByData(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_FilterByData()
        self.thisptr = <_Action_FilterByData*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_LESsplit(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_LESsplit()
        self.thisptr = <_Action_LESsplit*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_NativeContacts(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NativeContacts()
        self.thisptr = <_Action_NativeContacts*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_VelocityAutoCorr(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_VelocityAutoCorr()
        self.thisptr = <_Action_VelocityAutoCorr*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_SetVelocity(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_SetVelocity()
        self.thisptr = <_Action_SetVelocity*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_MultiVector(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MultiVector()
        self.thisptr = <_Action_MultiVector*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_MinImage(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MinImage()
        self.thisptr = <_Action_MinImage*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_ReplicateCell(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_ReplicateCell()
        self.thisptr = <_Action_ReplicateCell*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_AreaPerMol(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AreaPerMol()
        self.thisptr = <_Action_AreaPerMol*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Energy(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Energy()
        self.thisptr = <_Action_Energy*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_CheckChirality(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CheckChirality()
        self.thisptr = <_Action_CheckChirality*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Channel(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Channel()
        self.thisptr = <_Action_Channel*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()


cdef class Action_Volume(Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Volume()
        self.thisptr = <_Action_Volume*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()

