# distutils: language = c++
from cython.operator cimport dereference as deref

# TODO : double-check C++ code

cdef class ActionList:
    def __cinit__(self):
        self.thisptr = new _ActionList()

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def add_action(self, actionobj, 
                         ArgList arglist, 
                         top=None, 
                         FrameList flist=FrameList(), 
                         DataSetList dlist=DataSetList(), 
                         DataFileList dflist=DataFileList()):
        """
        Add action to ActionList

        Parameters:
        ==========
        actionobj :: Action object
        arglist :: ArgList instance
        toplist :: TopologyList instance
        flist :: FrameList instance
        dlist :: DataSetList 
        dflist :: DataFileList
        """
        cdef FunctPtr func = <FunctPtr> actionobj.alloc()
        cdef TopologyList toplist
        if isinstance(top, Topology):
            toplist = TopologyList()
            toplist.add_parm(top)
        elif isinstance(top, TopologyList):
            toplist = top
        # add function pointer: How?
        return self.thisptr.AddAction(func.ptr, arglist.thisptr[0], 
                                      toplist.thisptr, flist.thisptr, 
                                      dlist.thisptr, dflist.thisptr)

    def process(self, Topology top):
        # let cpptraj free mem
        top.py_free_mem = False
        return self.thisptr.SetupActions(&(top.thisptr))

    def do_actions(self, int idx=0, Frame frame=Frame()):
        # TODO : read cpptraj code to check memory stuff
        # set py_free_mem = False to let cpptraj does its job
        if frame.is_empty():
            raise ValueError("empty Frame, what can I do with this?")

        frame.py_free_mem = False
        return self.thisptr.DoActions(&(frame.thisptr), idx)

    def listinfo(self):
        self.thisptr.List()

    def is_empty(self):
        return self.thisptr.Empty()

    @property
    def n_actions(self):
        return self.thisptr.Naction()

    def cmd_string(self, int i):
        return self.thisptr.CmdString(i)

    def action_alloc(self, int i):
        # TODO : do we need to expose this method here?
        # return func_ptr
        cdef FunctPtr func = FunctPtr()
        if i >= self.n_actions:
            raise IndexError("index must be < " + str(self.n_actions)) 
        func.ptr = self.thisptr.ActionAlloc(i)
        return func
