from __future__ import absolute_import
from pytraj.AtomMask import AtomMask
from pytraj.externals.six import string_types
from pytraj.Frame import Frame


class FrameIter(object):
    """
    internal class
    create this class to hold all iterating information

    Parameters
    ----------
    top : new Topology
    original_top : original Topology
    start, stop, stride : int
    mask : str, defaul ""
        only take atom with given mask
    frame_indices : iterable, default: None
        if frame_indices is not None: ignore (start, stop, stride)
    autoimage : bool, default: False
        if autoimage, perform autoimage
    rmsfit : int or a tuple, defaul False
        if rmsfit, perform rms fit to reference. If ``rmsfit`` is an integer, perform
        rms fit to indicated frame for all atoms. If ``rmsfit`` is a tuple, perform rmsfit 
        to given frame with given mask. if bot ``autoimage`` and ``rmsfit`` are specified,
        do ``autoimage`` first
    n_frames : total number of frame. read-only
    copy : bool, defaul: True
        if True, always make a copy of Frame when iterating.

    Examples
    --------
    >>> # create FrameIter with start, stop, stride = 0, 8, 2
    >>> # autoimage=False, rmsfit=False
    >>> traj.iterframe(0, 8, 2)

    >>> # create FrameIter with start, stop, stride = 2, 8, 1
    >>> # autoimage=False, rmsfit=False
    >>> traj.iterframe(2, 8)

    >>> # create FrameIter with start, stop, stride = 2, 8, 1
    >>> # autoimage=False, rmsfit=False, mask='@CA'
    >>> traj.iterframe(2, 8, mask='@CA')

    >>> # create FrameIter with start, stop, stride = 2, 8, 1
    >>> # autoimage=True, rmsfit=False, mask='@CA'
    >>> traj.iterframe(2, 8, autoimage=True, mask='@CA')
    """

    def __init__(self, fi_generator,
                 original_top=None,
                 new_top=None,
                 start=0,
                 stop=-1,
                 stride=1,
                 mask="",
                 autoimage=False,
                 rmsfit=None,
                 is_trajiter=False,
                 n_frames=None,
                 copy=True,
                 frame_indices=None):
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
        self._n_frames = n_frames
        self.copy = copy
        self.frame_indices = frame_indices

    def __len__(self):
        return self._n_frames

    @property
    def n_frames(self):
        return self._n_frames

    @property
    def __name__(self):
        '''for inspecting
        '''
        return "FrameIter"

    def __str__(self):
        root_msg = '<FrameIter with '
        root_msg2 = 'start=%s, stop=%s, stride=%s, n_frames=%s, \n' % (
            self.start, self.stop, self.stride, self.n_frames)
        root_msg3 = 'frame_indices=%s \n' % self.frame_indices

        more_msg = 'autoimage=%s, rmsfit=%s, copy=%s> \n' % (
            self.autoimage, self.rmsfit, self.copy)
        return "".join((root_msg, root_msg2, root_msg3, more_msg))

    def __repr__(self):
        return self.__str__()

    def save(self,
             filename='',
             overwrite=False,
             mode='', *args, **kwd):
        '''save to different file format.

        Notes
        -----
        FrameIter will be exhausted since this is an iterator.

        Examples
        --------
        >>> fi = traj(2, 8, 2, mask='@CA')
        >>> fi.save('test.nc', overwrite=True)
        >>> # short version
        >>> traj(2, 8, 2, mask='@CA').save('test.nc')
        '''
        from pytraj.io import write_traj
        write_traj(filename=filename,
                   traj=self,
                   top=self.top,
                   indices=None,
                   overwrite=overwrite,
                   mode=mode, *args, **kwd)

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
            if self.is_trajiter and self.copy:
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
                        raise TypeError('')
                frame2 = Frame(atm.n_atoms)
                frame2.set_coords(frame, atm)

                yield frame2
            else:
                yield frame
