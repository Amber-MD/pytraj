from __future__ import absolute_import

#from .. Trajectory cimport Trajectory
#from .. Frame cimport _Frame, Frame
from .. Trajectory import Trajectory
from .. Frame import Frame
from .. externals import load_pseudo_parm
from .. _shared_methods import my_str_method

# TODO: should have Base class for all external packages (mdtraj, MDAnalysis, ...)
#cdef class TrajectoryMDAnalysisIterator(object):
class TrajectoryMDAnalysisIterator(object):
    #cdef object _traj_holder
    #cdef public Topology top

    #def __cinit__(self, mdanalysis_object, top=None):
    def __init__(self, mdanalysis_object, top=None):
        self._ext_holder = mdanalysis_object
        self._traj_holder = mdanalysis_object.trajectory

        if top is not None:
            self.top = top
        else:
            self.top = load_pseudo_parm(self._ext_holder)

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
       # cdef int n_atoms = self.n_atoms
        atom_groups = self._ext_holder.atoms

        for _ in self._traj_holder:
            frame = Frame(self.n_atoms)
            frame.xyz[:] = atom_groups.positions.astype('f8')
            yield frame

    def __getitem__(self, idx):
        i = 0
        atom_groups = self._ext_holder.atoms
        with self:
            for _ in self._traj_holder:
                if i == idx:
                    frame = Frame(self.n_atoms)
                    # for some reasons, need to set py_free_mem=Fale 
                    # to keep Frame's lifetime
                    frame.py_free_mem = False
                    frame.xyz[:] = atom_groups.positions
                    return frame # break the loop
                i += 1

    def __enter__(self):
        self._traj_holder.rewind()
        return self

    def __exit__(self, *args):
        self._traj_holder.close()

    @property
    def n_frames(self):
        return self._traj_holder.numframes

    @property
    def size(self):
        return self.n_frames

    @property
    def n_atoms(self):
        return self._traj_holder.numatoms
