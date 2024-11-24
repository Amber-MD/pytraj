# distutils: language = c++
from __future__ import print_function
from pytraj.utils.decorators import makesureABC
from pytraj.utils import is_generator
from pytraj.trajectory.shared_methods import iterframe_master
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
    ADICT["hbond"] = _ALL['Action_HydrogenBond']

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

        if isinstance(command, str):
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

        if status == ERR or status == SKIP:
            # cpptraj have a bunch of options, so we only check if there is
            # ERR or SKIP
            raise RuntimeError("Failed to setup action. Use pytraj._verbose() to "
                    "turn on the error report.")

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
        self.n_frames = 0 # Stop mark for generated script


cdef class BaseAction(Action):
    cdef _Action* baseptr
    cdef _Action* thisptr
    cdef bool own_memory

    def __cinit__(self, _Action* action_ptr):
        self.baseptr = <_Action*> action_ptr
        self.thisptr = <_Action*> self.baseptr
        self.own_memory = True

    def __dealloc__(self):
        if self.baseptr is not NULL and self.own_memory:
            del self.baseptr

    def help(self):
        self.thisptr.Help()

cdef class Action_Align(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Align())

cdef class Action_Angle(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Angle())

cdef class Action_AreaPerMol(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_AreaPerMol())

cdef class Action_AtomMap(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_AtomMap())

cdef class Action_AtomicCorr(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_AtomicCorr())

cdef class Action_AtomicFluct(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_AtomicFluct())

cdef class Action_AutoImage(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_AutoImage())

cdef class Action_Average(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Average())

cdef class Action_Bounds(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Bounds())

cdef class Action_Box(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Box())

cdef class Action_Center(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Center())

cdef class Action_Channel(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Channel())

cdef class Action_CheckChirality(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_CheckChirality())

cdef class Action_CheckStructure(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_CheckStructure())

cdef class Action_Closest(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Closest())

cdef class Action_ClusterDihedral(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_ClusterDihedral())

cdef class Action_Contacts(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Contacts())

cdef class Action_CreateCrd(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_CreateCrd())

cdef class Action_DNAionTracker(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_DNAionTracker())

cdef class Action_DSSP(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_DSSP())

cdef class Action_Density(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Density())

cdef class Action_Diffusion(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Diffusion())

cdef class Action_Dihedral(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Dihedral())

cdef class Action_Dipole(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Dipole())

cdef class Action_DistRmsd(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_DistRmsd())

cdef class Action_Distance(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Distance())

cdef class Action_Energy(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Energy())

cdef class Action_Esander(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Esander())

cdef class Action_FilterByData(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_FilterByData())

cdef class Action_FixAtomOrder(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_FixAtomOrder())

cdef class Action_FixImagedBonds(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_FixImagedBonds())

cdef class Action_GIST(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_GIST())

cdef class Action_Grid(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Grid())

cdef class Action_GridFreeEnergy(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_GridFreeEnergy())

cdef class Action_HydrogenBond(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_HydrogenBond())

cdef class Action_Image(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Image())

cdef class Action_Jcoupling(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Jcoupling())

cdef class Action_LESsplit(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_LESsplit())

cdef class Action_LIE(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_LIE())

cdef class Action_LipidOrder(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_LipidOrder())

cdef class Action_MakeStructure(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_MakeStructure())

cdef class Action_Mask(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Mask())

cdef class Action_Matrix(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Matrix())

cdef class Action_MinImage(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_MinImage())

cdef class Action_Molsurf(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Molsurf())

cdef class Action_MultiDihedral(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_MultiDihedral())

cdef class Action_MultiVector(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_MultiVector())

cdef class Action_NAstruct(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_NAstruct())

cdef class Action_NMRrst(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_NMRrst())

cdef class Action_NativeContacts(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_NativeContacts())

cdef class Action_OrderParameter(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_OrderParameter())

cdef class Action_Outtraj(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Outtraj())

cdef class Action_PairDist(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_PairDist())

cdef class Action_Pairwise(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Pairwise())

cdef class Action_Principal(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Principal())

cdef class Action_Projection(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Projection())

cdef class Action_Pucker(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Pucker())

cdef class Action_Radgyr(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Radgyr())

cdef class Action_Radial(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Radial())

cdef class Action_RandomizeIons(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_RandomizeIons())

cdef class Action_Remap(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Remap())

cdef class Action_ReplicateCell(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_ReplicateCell())

cdef class Action_Rmsd(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Rmsd())

cdef class Action_Rotate(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Rotate())

cdef class Action_RunningAvg(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_RunningAvg())

cdef class Action_STFC_Diffusion(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_STFC_Diffusion())

cdef class Action_Scale(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Scale())

cdef class Action_SetVelocity(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_SetVelocity())

cdef class Action_Spam(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Spam())

cdef class Action_Strip(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Strip())

cdef class Action_Surf(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Surf())

cdef class Action_SymmetricRmsd(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_SymmetricRmsd())

cdef class Action_Temperature(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Temperature())

cdef class Action_Translate(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Translate())

cdef class Action_Unstrip(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Unstrip())

cdef class Action_Unwrap(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Unwrap())

cdef class Action_Vector(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Vector())

cdef class Action_VelocityAutoCorr(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_VelocityAutoCorr())

cdef class Action_Volmap(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Volmap())

cdef class Action_Volume(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Volume())

cdef class Action_Watershell(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_Watershell())

cdef class Action_XtalSymm(BaseAction):
    def __cinit__(self):
        BaseAction.__cinit__(self, <_Action*> new _Action_XtalSymm())
