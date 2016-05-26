from nglview import NGLWidget, PyTrajTrajectory

class TrajectoryViewer(NGLWidget):

    def add_pytraj_trajectory(self, traj):
        """add pytraj's Trajectory or TrajectoryIterator object
        """
        self.add_trajectory(PyTrajTrajectory(traj))
