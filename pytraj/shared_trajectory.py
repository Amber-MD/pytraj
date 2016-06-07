from .shared_methods import my_str_method

class SharedTrajectory(object):

    def view(self, *args, **kwargs):
        """require NGLView

        Parameters
        ----------
        args and kwargs : NGLView's arguments
        """
        from pytraj.visualization.nglview_wrapper import TrajectoryViewer
        from nglview import PyTrajTrajectory

        return TrajectoryViewer(PyTrajTrajectory(self), *args, **kwargs)

    def save(self,
             filename="",
             overwrite=False,
             **kwd):
        from pytraj.c_traj.c_trajout import TrajectoryWriter
    
        with TrajectoryWriter(filename=filename,
                     top=self.top,
                     overwrite=overwrite,
                     **kwd) as trajout:

            for frame in self:
                trajout.write(frame)

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return str(self)
