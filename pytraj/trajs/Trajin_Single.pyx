# distutils: language = c++
import os

cdef class Trajin_Single(Trajin):
    def __cinit__(self, filename=None, top=None, *args, **kwd):
        self.baseptr_1 = <_Trajin*> new _Trajin_Single()
        self.thisptr = <_Trajin_Single*> self.baseptr_1

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
        filename = filename.encode("UTF-8")
        if not top.is_empty():
            #print "update Topology for %s instance" % (self.__class__.__name__)
            self._top = top.copy()
        return self.thisptr.SetupTrajRead(filename, arglist.thisptr[0], self._top.thisptr, check_box)
