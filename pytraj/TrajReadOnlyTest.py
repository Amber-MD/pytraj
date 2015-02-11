# distutils: language = c++
from __future__ import print_function
import os
from pytraj.externals.six import string_types
from pytraj.Trajin_Single import Trajin_Single
from pytraj.TrajinList import TrajinList

class TrajReadOnlyTest():
    def __init__(self, filename=None, top=None, *args):
        self._trajlist = []
        self.top = top
        if filename:
            if top:
                if isinstance(top, string_types):
                    top_ = Topology(top)
                elif isinstance(top, Topology):
                    top_ = top.copy()
                if not self.top.is_empty():
                    print ("Repalce self.top with new provided top")
                self.top = top_.copy()
            self.load(filename)
        
    def __iter__(self):
        for traj in self._trajlist:
            for frame in traj:
                yield frame

    def copy(self):
        other = TrajReadOnlyTest()
        other.top = self.top.copy()
        other._trajlist = self._trajlist[:]
        return other

    # Let base-class Trajin take care those methods?
    def load(self, *args, **kwd):
        """
        Load trajectory from file.

        Parameters:
        filename :: string (trajectory file's name)
        ArgList instance
        Topology instance
        chexbox :: (default = True)
        """
        _traj = Trajin_Single()
        _traj.top = self.top.copy()
        # Currently we can not assigne self.top to top.copy() since Cython does not know self.top type
        # need to use self._top since we declare it in TrajectoryFile.pxd
        _traj.load(*args, **kwd)
        self._trajlist.append(_traj)

    @property
    def n_frames(self):
        n = 0
        for traj in self._trajlist:
            n += traj.n_frames
        return n

    @property
    def n_trajs(self):
        return len(self._trajlist)

    def join(self, other):
        if self == other:
            raise ValueError("join yourself, noway")
        if isinstance(other, Trajin_Single):
            self._trajlist.append(other)
        elif isinstance(other, TrajReadOnlyTest):
            for traj in other._trajlist:
                self._trajlist.append(traj)
        elif isinstance(other, (TrajinList, tuple, list)):
            # list of Trajin
            for traj in other:
                self._trajlist.append(traj)
            
    def is_empty(self):
        return self.n_frames == 0

    def strip_atoms(self, *args, **kwd):
        raise NotImplementedError("Read-only-traj, use FrameArray class")
