# distutils: language = c++

# this code was translated from cpptraj (written in C++ by Daniel R. Roe)
# Translation: Hai Nguyen

from copy import copy
from libcpp.string cimport string
from ..Box cimport Box
from ..AtomMask cimport AtomMask
from ..Topology cimport Topology
from ..Frame cimport Frame
from ..ParmFile cimport *

cdef class Action_Strip:
    cdef public Topology old_parm, new_parm
    cdef public Frame new_frame
    cdef public bint nobox  
    cdef public string prefix
    cdef public AtomMask M1

    def __cinit__(self, nobox=False, prefix="test"):
        self.old_parm = Topology()
        self.new_parm = Topology()
        self.new_frame = Frame()
        self.nobox = nobox
        self.prefix = prefix 
        self.M1 = AtomMask()
        
    def init(self,ArgList arglist, debug=0):
        self.prefix = arglist.get_string_key("outprefix")
        self.nobox = arglist.has_key("nobox")
        mask1 = arglist.get_mask_next()
        self.M1.set_mask_string(mask1)
        self.M1.invert_mask()

    def setup(self, Topology current_parm):
        if current_parm.setup_integer_mask(self.M1):
            return "Error"
        s_num = current_parm.n_atoms - self.M1.n_selected
        self.old_parm.thisptr[0] = current_parm.thisptr[0]
        try:
            self.new_parm = current_parm.modify_state_by_mask(self.M1)
        except:
            raise ValueError()

        if self.nobox:
            self.new_parm.set_box(Box())

        self.new_frame.set_frame_v(self.new_parm, self.new_parm.has_vel, self.new_parm.n_repdim)

    def do_action(self, frame_num, Frame current_frame):
        self.new_frame.set_frame(current_frame, self.M1)

    def write_parm(self, top=None, prefix="temp", parmtype=AMBERPARM):
        cdef ParmFile pfile = ParmFile()
        pfile.write_prefix_topology(<Topology>top, <string>prefix, parmtype, 0)
