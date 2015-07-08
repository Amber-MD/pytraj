from __future__ import absolute_import
from pytraj.AtomMask import AtomMask
from pytraj.compat import string_types
from pytraj.Frame import Frame


class FrameIter(object):

    """
    internal class, not for user
    create this class to hold topology
    only iterate this once
    """

    def __init__(self, fi_generator,
                 original_top=None,
                 new_top=None, start=0, stop=-1, stride=1,
                 mask="", autoimage=False,
                 rmsfit=None,
                 is_trajiter=False,
                 n_frames=None):
        self.top = new_top
        self.original_top = original_top
        self.frame_iter = fi_generator
        self.start = start
        self.stop = stop
        self.stride = stride
        self.mask = mask
        self.autoimage = autoimage
        self.rmsfit = rmsfit
        # use `copy_frame` for TrajectoryIterator
        self.is_trajiter = is_trajiter
        self.n_frames = n_frames

    @property
    def __name__(self):
        return "FrameIter"

    def __str__(self):
        root_msg = '<pytraj.core.frameiter.FrameIter with '
        root_msg2 = 'start=%s, stop=%s, stride=%s \n' % (
            self.start, self.stop, self.stride)

        more_msg = 'autoimage=%s, rmsfit=%s> \n' % (
            self.autoimage, self.rmsfit)
        return "".join((root_msg, root_msg2, more_msg))

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        if self.autoimage:
            from pytraj.actions.CpptrajActions import Action_AutoImage
            image_act = Action_AutoImage()
            image_act.read_input("", top=self.original_top)
            image_act.process(self.original_top)
        if self.rmsfit is not None:
            try:
                ref, mask_for_rmsfit = self.rmsfit
            except ValueError:
                ref = self.rmsfit[0]
                mask_for_rmsfit = "*"
            need_align = True
            if self.autoimage:
                # need to do autoimage for ref too
                # make a copy to avoid changing ref
                ref = ref.copy()
                image_act.do_action(ref)
            from pytraj.actions.CpptrajActions import Action_Rmsd
            rmsd_act = Action_Rmsd()
            rmsd_act.read_input(mask_for_rmsfit, top=self.original_top)
            rmsd_act.process(self.original_top)
            # creat first frame to trick cpptraj to align to this.
            rmsd_act.do_action(ref)
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for frame0 in self.frame_iter:
            if self.is_trajiter:
                # use copy for TrajectoryIterator
                # so [f for f in traj()] will return a list of different 
                # frames
                frame = frame0.copy()
            else:
                frame = frame0
            if self.autoimage:
                #from pytraj.actions.CpptrajActions import Action_AutoImage
                #Action_AutoImage()("", frame, self.top)
                image_act.do_action(frame)
            if need_align:
                # trick cpptraj to fit to 1st frame (=ref)
                rmsd_act.do_action(frame)
            if self.mask is not None:
                mask = self.mask
                if isinstance(mask, string_types):
                    atm = self.original_top(mask)
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
