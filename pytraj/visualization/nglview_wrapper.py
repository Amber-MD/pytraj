from nglview import NGLWidget, PyTrajTrajectory

class TrajectoryViewer(NGLWidget):

    def add(self, traj, *args, **kwargs):
        """add pytraj's Trajectory or TrajectoryIterator object
        """
        from pytraj import Trajectory, TrajectoryIterator

        if isinstance(traj, (Trajectory, TrajectoryIterator)):
            super(TrajectoryViewer, self).add_trajectory(PyTrajTrajectory(traj), *args, **kwargs)
        else:
            # try with your own risk
            super(TrajectoryViewer, self).add_trajectory(traj, *args, **kwargs)

    def add_trajectory(self, traj, *args, **kwargs):
        self.add(traj, *args, **kwargs)
