from nglview import NGLWidget, PyTrajTrajectory

class TrajectoryViewer(NGLWidget):

    def add(self, traj, *args, **kwargs):
        """add pytraj's Trajectory or TrajectoryIterator object
        """
        from pytraj import Trajectory, TrajectoryIterator

        if isinstance(traj, (Trajectory, TrajectoryIterator)):
            self.add_trajectory(PyTrajTrajectory(traj), *args, **kwargs)
        else:
            # try with your own risk
            self.add_component(traj, *args, **kwargs)
