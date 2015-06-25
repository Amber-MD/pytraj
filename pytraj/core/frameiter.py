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
                 rmsfit=None):
        self.new_top = new_top
        self.original_top = original_top
        self.frame_iter = fi_generator
        self.start = start
        self.stop = stop
        self.stride = stride
        self.mask = mask
        self.autoimage = autoimage
        self.rmsfit = rmsfit

    @property
    def __name__(self):
        return "FrameIter"

    def __str__(self):
        root_msg = '<pytraj.core.frameiter.FrameIter with '
        root_msg2 = 'start=%s, stop=%s, stride=%s \n' % (self.start, self.stop, self.stride)

        more_msg = 'autoimage=%s, rmsfit=%s> \n' % (self.autoimage, self.rmsfit) 
        return "".join((root_msg, root_msg2, more_msg))

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        if self.autoimage:
            from pytraj.actions.CpptrajActions import Action_AutoImage
            act = Action_AutoImage()
        if self.rmsfit is not None:
            try:
                ref, mask_for_rmsfit = self.rmsfit
            except ValueError:
                ref = self.rmsfit[0]
                mask_for_rmsfit = "*"
            need_align = True
            from pytraj.actions.CpptrajActions import Action_Rmsd
        else:
            need_align = False
            ref, mask_for_rmsfit = None, None

        for frame in self.frame_iter:
            if self.autoimage:
                act(current_frame=frame, top=self.original_top)
            if need_align:
                act = Action_Rmsd()
                # trick cpptraj to fit to 1st frame (=ref)
                act(mask_for_rmsfit, [ref, frame], top=self.original_top)
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
