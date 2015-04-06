import os
import unittest
from pytraj.trajs.TrajectoryFile import TrajectoryFile

class TestTraj(unittest.TestCase):
    def test_1(self):
        traj = TrajectoryFile()
        t = traj._get_type_from_ext(".nc")
        assert t == "AMBERNETCDF"
        print(t)
        ext = traj._get_ext_for_type("AMBERNETCDF")
        assert ext == ".nc"
        print(ext)
        fs = traj._format_string("AMBERNETCDF")
        assert fs == "Amber NetCDF"
        print(fs)

    def test_classmethod(self):
        t = TrajectoryFile._get_type_from_ext(".nc")
        assert t == "AMBERNETCDF"
        print(t)
        ext = TrajectoryFile._get_ext_for_type("AMBERNETCDF")
        assert ext == ".nc"
        print(ext)
        fs = TrajectoryFile._format_string("AMBERNETCDF")
        assert fs == "Amber NetCDF"
        print(fs)
        print(dir(TrajectoryFile))
        TrajectoryFile._read_options()
        TrajectoryFile._write_options()

if __name__ == "__main__":
    unittest.main()
