# distutils: language = c++
import os

cdef class Trajin_Single(Trajin):
    def __cinit__(self, filename=None, top=None, *args, **kwd):
        # thisptr is from Trajin class
        # now it points to derived class
        self.baseptr0 = <_TrajectoryFile*> new _Trajin_Single()
        self.baseptr_1 = <_Trajin*> self.baseptr0
        self.thisptr = <_Trajin_Single*> self.baseptr0

        if filename:
            if top:
                if isinstance(top, basestring):
                    top_ = Topology(top)
                elif isinstance(top, Topology):
                    top_ = top.copy()
                if not self.top.is_empty():
                    print "Repalce self.top with new provided top"
                self.top = top_.copy()
            self.load(filename, *args, **kwd)
        
    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def alloc(self):
        # TODO : get Segmentation fault error when:
        # for traj in trajin:
        #    pass
        """return Trajin view"""
        cdef Trajin trajin = Trajin()
        trajin.baseptr_1 = <_Trajin*> self.thisptr
        # re-cast baseptr0 too
        # anyway to avoid re_casting?
        trajin.baseptr0 = <_TrajectoryFile*> self.thisptr
        # since we just get a view of `self`, we let `self` free memory for self.top
        # don't let trajin.top do this
        return trajin

    # Let base-class Trajin take care those methods?
    def load(Trajin_Single self, filename='', Topology top=Topology(), 
             ArgList arglist=ArgList(), bint check_box=True,
             ):
        """
        Load trajectory from file.

        Parameters:
        filename :: string (trajectory file's name)
        ArgList instance
        Topology instance
        chexbox :: (default = True)
        """
        # Currently we can not assigne self.top to top.copy() since Cython does not know self.top type
        # need to use self._top since we declare it in TrajectoryFile.pxd
        filename = filename.encode("UTF-8")
        if not top.is_empty():
            #print "update Topology for %s instance" % (self.__class__.__name__)
            self._top = top.copy()
        return self.thisptr.SetupTrajRead(filename, arglist.thisptr[0], self._top.thisptr, check_box)
