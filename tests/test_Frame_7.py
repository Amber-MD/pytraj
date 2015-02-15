import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj.externals import netcdf
        top = Topology("./data/tz2.ortho.parm7")
        fh = netcdf.netcdf_file("./data/tz2.ortho.nc")
        print (fh)
        traj = mdio.load("./data/tz2.ortho.nc", top)
        coords = fh.variables['coordinates']
        farr0 = FrameArray()
        farr0.top = top.copy()
        for arr0 in coords:
            frame = Frame(top.n_atoms)
            frame.set_from_crd(arr0.flatten())
            farr0.append(frame)
        print (farr0.size)

        for i in range(traj.size):
            assert traj[i].same_coords_as(farr0[i]) == True
        fh.close()
            
if __name__ == "__main__":
    unittest.main()
