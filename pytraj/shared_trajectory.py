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
