from __future__ import absolute_import

#from .. Trajectory cimport Trajectory
#from .. Frame cimport _Frame, Frame
from .. Trajectory import Trajectory
from .. Frame import Frame
from .. externals import load_pseudo_parm

#cdef class TrajectoryMDAnalysisIterator(object):
class TrajectoryMDAnalysisIterator(object):
    #cdef object _traj_holder
    #cdef public Topology top

    #def __cinit__(self, mdanalysis_object, top=None):
    def __init__(self, mdanalysis_object, top=None):
        self._traj_holder = mdanalysis_object.trajectory
        if top is not None:
            self.top = top
        else:
            self.top = load_pseudo_parm(mdanalysis_object)

    def __iter__(self):
       # cdef int n_atoms = self.n_atoms
        atom_groups = self._traj_holder.atoms

        for _ in self._traj_holder.trajectory:
            frame = Frame(self.n_atoms)
            frame.xyz[:] = atom_groups.positions.astype('f8')
            yield frame

    @property
    def n_frames(self):
        return self._traj_holder.numframes

    @property
    def n_atoms(self):
        return self._traj_holder.numatoms
