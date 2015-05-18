from __future__ import absolute_import

#from .. Trajectory cimport Trajectory
#from .. Frame cimport _Frame, Frame
from .. utils.check_and_assert import is_int
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
        atom_groups = self._ext_holder.atoms
        with self:
            if is_int(idx):
                if idx >= self.n_frames or idx < 0: 
                    raise ValueError("must have 0 <= idx < self.n_frames")
                i = 0
                for _ in self._traj_holder:
                    if i == idx:
                        frame = Frame(self.n_atoms)
                        # for some reasons, need to set py_free_mem=Fale 
                        # to keep Frame's lifetime
                        frame.py_free_mem = False
                        frame.xyz[:] = atom_groups.positions
                        return frame # break the loop
                    i += 1
            elif isinstance(idx, slice):
                fa = Trajectory()
                fa.top = self.top
                atom_groups = self._ext_holder.atoms
                start, stop, stride = idx.indices(self.n_frames)

                try:
                    # if MDAnalysis object support slicing
                    for _ in self._traj_holder[idx]:
                        frame = Frame(self.n_atoms)
                        frame.xyz[:] = atom_groups.positions
                        fa.append(frame, copy=False)
                    return fa
                except:
                    # old fashion way
                    count = start
                    while count < stop:
                        fa.append(self[count], copy=False)
                        count += stride
                    return fa
            else:
                raise NotImplementedError()

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

    @property
    def xyz(self):
        from pytraj.io import load_MDAnalysis
        return load_MDAnalysis(self._ext_holder, top=self.top).xyz
