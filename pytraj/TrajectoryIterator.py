"""This is a thin wrapper of Trajin_Single
We need to sub-class Trajin_Single to use Trajectory
(we called Trajin_Single from Trajectory, so we can not call Trajectory back from 
Trajin_Single)
"""
from __future__ import absolute_import
from pytraj.utils.check_and_assert import is_pytraj_trajectory
from pytraj.trajs.Trajin_Single import Trajin_Single
from pytraj._action_in_traj import ActionInTraj
from pytraj.action_dict import ActionDict
from pytraj.Frame import Frame
from pytraj.AtomMask import AtomMask
from pytraj.externals.six import string_types
from pytraj.exceptions import PytrajMemviewError


class TrajectoryIterator(Trajin_Single, ActionInTraj):
    def __init__(self, *args, **kwd):
        pass

    @property
    def topology(self):
        """traditional name for Topology file"""
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop

    def frame_iter(self, start=0, stop=-1, stride=1, mask=None, autoimage=False, rmsfit_to=None):
        """frame iterator

        Parameters
        ----------
        start : int, default=0
        stop : int, default=-1 (last frame)
        stride : int, deault=1 (no skipping)
        mask : str | array of integers | AtomMask object, default=None
            get sub-frame with coords of given mask
        autoimage : bool, default=False
            perform `autoimage`
        rmsfit_to : None | tuple/list of (Frame, mask)
            perform `rmsfit_to` reference if not `None`
            Notes : reference must have the same n_atoms
        """
        if autoimage:
            act = ActionDict()['autoimage']
        if rmsfit_to is not None:
            ref, mask_for_rmsfit = rmsfit_to
            need_align = True
            from pytraj.actions.Action_Rmsd import Action_Rmsd
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for frame in super(TrajectoryIterator, self).frame_iter(start, stop, stride):
            if autoimage:
                act(current_frame=frame, top=self.top)
            if need_align:
                act = Action_Rmsd()
                # trick cpptraj to fit to 1st frame (=ref)
                act(mask_for_rmsfit, [ref, frame], top=self.top)
            if mask is not None:
                if isinstance(mask, string_types):
                    atm = self.top(mask)
                else:
                    try:
                        atm = AtomMask()
                        atm.add_selected_indices(mask)
                    except TypeError:
                        raise PytrajMemviewError()
                frame2 = Frame(atm.n_atoms)
                frame2.set_coords(frame, atm)
                yield frame2
            else:
                yield frame

    def chunk_iter(self, chunk=2, start=0, stop=-1, 
                   autoimage=False, 
                   rmsfit_to=None,
                   copy_top=False):
        """
        Parameters
        ----------
        chunk : int, default=2
            size of each chunk. Notes: final chunk's size might be changed
        start : int, default=0 (first frame)
        stop : int, default=-1 (last frame)
        autoimage : bool, default=False
        rmsfit_to : None | tuple/list of (reference frame, mask)

        Examples
        --------
            for chunk in trajiter.chunk_iter(100, autoimage=True, rmsfit_to=(ref0, '@CA'))
        """

        if rmsfit_to is not None:
            ref, mask_for_rmsfit = rmsfit_to
            need_align = True
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for chunk in super(TrajectoryIterator, self).chunk_iter(chunk, start, stop, copy_top):
            # always perform autoimage before doing fitting
            # chunk is `Trajectory` object, having very fast `autoimage` and `rmsfit_to` methods
            if autoimage:
                chunk.autoimage()
            if need_align:
                chunk.rmsfit_to(ref, mask_for_rmsfit)
            yield chunk
