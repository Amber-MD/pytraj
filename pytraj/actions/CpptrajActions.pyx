# distutils: language = c++
from __future__ import print_function
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types
from pytraj.utils import is_generator
from pytraj.utils.check_and_assert import is_pytraj_trajectory
from pytraj._shared_methods import iterframe_master
from cython.operator cimport dereference as deref


cdef class Action:
    '''interface to Cpptraj's Action
    '''
    def __cinit__(self):
        # don't directly create instance of this ABC class.
        self.n_frames = 0
        self.top_is_processed = False

    def __dealloc__(self):
        # should I del pointer here or in subclass? 
        #del self.baseptr
        pass

    def __del__(self):
        del self.baseptr

    def __str__(self):
        txt = "< %s object >" % (self.__class__.__name__)
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
                   top=TopologyList(),
                   DatasetList dslist=DatasetList(), 
                   DataFileList dflist=DataFileList(), 
                   int debug=0):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Action_Box().help())
        top : Topology or TopologyList instance, default=TopologyList()
        #flist : FrameList instance, default=FrameList()
        dslist : DatasetList instance, default=DatasetList()
        dflist : DataFileList instance, default=DataFileList()
        debug : int, default=0
            debug option from cpptraj. (Do we need this?)
        """
        cdef ArgList arglist
        cdef TopologyList toplist
        cdef RetType i_fail

        if isinstance(top, Topology):
            toplist = TopologyList()
            toplist.add_parm(top)
        elif isinstance(top, TopologyList):
            toplist = <TopologyList> top
        self.toplist = toplist

        if isinstance(command, string_types):
            #command = command.encode("UTF-8")
            arglist = ArgList(command)
        elif isinstance(command, ArgList):
            arglist = <ArgList> command

        i_fail = self.baseptr.Init(arglist.thisptr[0], toplist.thisptr, 
                       dslist.thisptr, dflist.thisptr,
                       debug)

        if i_fail != OK:
            # check before do_action to avoid segfault
            raise ValueError("")
        else:
            return i_fail

    @makesureABC("Action")
    def process(self, Topology top=Topology(), Topology new_top=Topology()): 
        """
        Process input and do initial setup
        (TODO : add more doc)

        Parameters:
        ----------
        top : Topology instance, default (no default)
        new_top : new Topology instance, default=Topology()
            Need to provide this instance if you want to change topology
        """
        if "Strip" in self.__class__.__name__:
            # since `Action_Strip` will copy a modified version of `top` and 
            # store in new_top, then __dealloc__ (from cpptraj)
            # we need to see _own_memory to False
            new_top._own_memory = False
        self.top_is_processed = True
        return self.baseptr.Setup(top.thisptr, &(new_top.thisptr))

    @makesureABC("Action")
    def do_action(self, current_frame=None, Frame new_frame=Frame(), update_mass=True):
        """
        Perform action on Frame. 
        Parameters:
        ----------
        current_frame : Frame instance need to be processed, default=Frame() 

        new_frame : Frame instance, defaul=Frame()
            if action change Frame, you need to have this
        """
        # debug
        cdef Frame frame
        cdef int i
        cdef object traj, tmptraj, farray

        if self.__class__.__name__ == 'Action_Strip':
            # let cpptraj do its job for this special action
            new_frame._own_memory = False

        if isinstance(current_frame, Frame):
            frame = <Frame> current_frame
            # make sure to update frame mass
            if update_mass:
                frame.set_frame_mass(self.toplist[0])
            self.baseptr.DoAction(self.n_frames, frame.thisptr, &(new_frame.thisptr))
            self.n_frames += 1
        else:
            for frame in iterframe_master(current_frame):
                self.do_action(frame, new_frame)

    @makesureABC("Action")
    def print_output(self):
        """Do we need this?"""
        self.baseptr.Print()

    # Do we really need this method?
    @classmethod
    def _get_action_from_functptr(cls, FunctPtr funct):
        cdef Action act = Action()
        if funct.ptr() == NULL:
            raise ValueError("NULL pointer")
        act.baseptr = <_Action*> funct.ptr()
        return act

    def _master(self, command='',
                  current_frame=Frame(),
                  top=Topology(),
                  dslist=DatasetList(), 
                  dflist=DataFileList(), 
                  new_top=Topology(),
                  new_frame=Frame(),
                  update_mass=True,
                  int debug=0,
                  update_frame=False,
                  quick_get=False):
        """
        TODO : (do we need this method?)
            + add doc
            + don't work with `chunk_iter`

        """
        if top.is_empty():
            _top = current_frame.top
        else:
            _top = top
        self.read_input(command=command, 
                        top=_top, 
                        dslist=dslist,
                        dflist=dflist, debug=debug)

        self.process(top=_top, new_top=new_top)
        self.do_action(current_frame, new_frame, update_mass=update_mass)
        return dslist

    def reset_counter(self):
        self.n_frames = 0


cdef class Action_Angle (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Angle()
        self.thisptr = <_Action_Angle*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_AreaPerMol (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AreaPerMol()
        self.thisptr = <_Action_AreaPerMol*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_AtomMap (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomMap()
        self.thisptr = <_Action_AtomMap*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_AtomicCorr (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomicCorr()
        self.thisptr = <_Action_AtomicCorr*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_AtomicFluct (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AtomicFluct()
        self.thisptr = <_Action_AtomicFluct*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_AutoImage (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_AutoImage()
        self.thisptr = <_Action_AutoImage*> self.baseptr

    def __dealloc__(self):
        if self.baseptr:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Average (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Average()
        self.thisptr = <_Action_Average*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Bounds (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Bounds()
        self.thisptr = <_Action_Bounds*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Box (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Box()
        self.thisptr = <_Action_Box*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Center (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Center()
        self.thisptr = <_Action_Center*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Channel (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Channel()
        self.thisptr = <_Action_Channel*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_CheckChirality (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CheckChirality()
        self.thisptr = <_Action_CheckChirality*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_CheckStructure (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CheckStructure()
        self.thisptr = <_Action_CheckStructure*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Closest (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Closest()
        self.thisptr = <_Action_Closest*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_ClusterDihedral (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_ClusterDihedral()
        self.thisptr = <_Action_ClusterDihedral*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Contacts (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Contacts()
        self.thisptr = <_Action_Contacts*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_CreateCrd (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_CreateCrd()
        self.thisptr = <_Action_CreateCrd*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_DNAionTracker (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DNAionTracker()
        self.thisptr = <_Action_DNAionTracker*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_DSSP (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DSSP()
        self.thisptr = <_Action_DSSP*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Density (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Density()
        self.thisptr = <_Action_Density*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Diffusion (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Diffusion()
        self.thisptr = <_Action_Diffusion*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()

cdef class Action_Dihedral (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Dihedral()
        self.thisptr = <_Action_Dihedral*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func

    def help(self):
        self.thisptr.Help()


cdef class Action_DihedralScan (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DihedralScan()
        self.thisptr = <_Action_DihedralScan*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Dipole (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Dipole()
        self.thisptr = <_Action_Dipole*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_DistRmsd (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_DistRmsd()
        self.thisptr = <_Action_DistRmsd*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Distance (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Distance()
        self.thisptr = <_Action_Distance*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()

    #def NOE_args(self, ArgList argIn, double noe_lbound, double noe_ubound, double noe_rexp):
    #    return self.thisptr.NOE_Args(argIn.thisptr[0], noe_lbound, noe_ubound, noe_rexp)



cdef class Action_Energy (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Energy()
        self.thisptr = <_Action_Energy*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_FilterByData (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_FilterByData()
        self.thisptr = <_Action_FilterByData*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_FixAtomOrder (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_FixAtomOrder()
        self.thisptr = <_Action_FixAtomOrder*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Gist (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Gist()
        self.thisptr = <_Action_Gist*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Grid (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Grid()
        self.thisptr = <_Action_Grid*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_GridFreeEnergy (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_GridFreeEnergy()
        self.thisptr = <_Action_GridFreeEnergy*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Hbond (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Hbond()
        self.thisptr = <_Action_Hbond*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Image (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Image()
        self.thisptr = <_Action_Image*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Jcoupling (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Jcoupling()
        self.thisptr = <_Action_Jcoupling*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_LESsplit (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_LESsplit()
        self.thisptr = <_Action_LESsplit*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_LIE (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_LIE()
        self.thisptr = <_Action_LIE*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_MakeStructure (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MakeStructure()
        self.thisptr = <_Action_MakeStructure*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Mask (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Mask()
        self.thisptr = <_Action_Mask*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Matrix (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Matrix()
        self.thisptr = <_Action_Matrix*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_MinImage (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MinImage()
        self.thisptr = <_Action_MinImage*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Molsurf (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Molsurf()
        self.thisptr = <_Action_Molsurf*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_MultiDihedral (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MultiDihedral()
        self.thisptr = <_Action_MultiDihedral*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_MultiVector (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_MultiVector()
        self.thisptr = <_Action_MultiVector*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_NAstruct (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NAstruct()
        self.thisptr = <_Action_NAstruct*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_NMRrst (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NMRrst()
        self.thisptr = <_Action_NMRrst*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_NativeContacts (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_NativeContacts()
        self.thisptr = <_Action_NativeContacts*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_OrderParameter (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_OrderParameter()
        self.thisptr = <_Action_OrderParameter*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Outtraj (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Outtraj()
        self.thisptr = <_Action_Outtraj*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_PairDist (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_PairDist()
        self.thisptr = <_Action_PairDist*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Pairwise (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Pairwise()
        self.thisptr = <_Action_Pairwise*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Principal (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Principal()
        self.thisptr = <_Action_Principal*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Projection (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Projection()
        self.thisptr = <_Action_Projection*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Pucker (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Pucker()
        self.thisptr = <_Action_Pucker*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Radgyr (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Radgyr()
        self.thisptr = <_Action_Radgyr*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Radial (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Radial()
        self.thisptr = <_Action_Radial*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_RandomizeIons (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_RandomizeIons()
        self.thisptr = <_Action_RandomizeIons*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_ReplicateCell (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_ReplicateCell()
        self.thisptr = <_Action_ReplicateCell*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Rmsd (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Rmsd()
        self.thisptr = <_Action_Rmsd*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func

    def help(self):
        self.thisptr.Help()


cdef class Action_Rotate (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Rotate()
        self.thisptr = <_Action_Rotate*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_RunningAvg (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_RunningAvg()
        self.thisptr = <_Action_RunningAvg*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_STFC_Diffusion (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_STFC_Diffusion()
        self.thisptr = <_Action_STFC_Diffusion*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Scale (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Scale()
        self.thisptr = <_Action_Scale*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_SetVelocity (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_SetVelocity()
        self.thisptr = <_Action_SetVelocity*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Spam (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Spam()
        self.thisptr = <_Action_Spam*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Strip (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Strip()
        self.thisptr = <_Action_Strip*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Surf (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Surf()
        self.thisptr = <_Action_Surf*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_SymmetricRmsd (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_SymmetricRmsd()
        self.thisptr = <_Action_SymmetricRmsd*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Temperature (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Temperature()
        self.thisptr = <_Action_Temperature*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Translate (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Translate()
        self.thisptr = <_Action_Translate*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Unwrap (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Unwrap()
        self.thisptr = <_Action_Unwrap*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Vector (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Vector()
        self.thisptr = <_Action_Vector*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_VelocityAutoCorr (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_VelocityAutoCorr()
        self.thisptr = <_Action_VelocityAutoCorr*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Volmap (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Volmap()
        self.thisptr = <_Action_Volmap*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Volume (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Volume()
        self.thisptr = <_Action_Volume*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()


cdef class Action_Watershell (Action):
    def __cinit__(self):
        self.baseptr = <_Action*> new _Action_Watershell()
        self.thisptr = <_Action_Watershell*> self.baseptr

    def __dealloc__(self):
        if self.baseptr is not NULL:
            del self.baseptr

    def alloc(self):
        """return a function-pointer object to be used with ActionList class
        """
        cdef FunctPtr func = FunctPtr()
        func.ptr = &(self.thisptr.Alloc)
        return func
        
    def help(self):
        self.thisptr.Help()
