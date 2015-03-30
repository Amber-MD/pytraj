import unittest
from pytraj.io import load
from pytraj import adict

class TestRandomizeIons(unittest.TestCase):
    def test_0(self):
        # get `traj` instance (FrameArray)
        traj = load(filename="./Test_RandomizeIons/adh206.tip3p.rst7.gz", 
                    top="./Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")
        top = traj.top.copy()
        # get 1st frame from `traj`
        frame0 = traj[0]
        frame1 = frame0.copy()
        bkframe = traj[0].copy()
        # make sure to reproduce cpptraj output
        savedframe = load(filename="./Test_RandomizeIons/random.crd.save",
                          top="./Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")[0]

        act1 = adict['randomizeions']
        act1("randomizeions @Na+ around :1-16 by 5.0 overlap 3.0", frame0, traj.top)
        assert frame0.rmsd(savedframe) < 1E-3
        # make sure changing frame0 won't change frame1
        assert frame0.rmsd(frame1) > 1E-1
        print ('rmsd frame0 and frame0 = ', frame0.rmsd(frame1))

        # try randomizeions for frame1
        # FIXME : can not continuously use more than 1 action
        from pytraj import allactions
        act2 = allactions.Action_RandomizeIons()
        print (hex(id(act1)))
        print (hex(id(act2)))

        act2("randomizeions @Na+ around :1-16 by 5.0 overlap 3.0", frame1, top)
        print ('rmsd frame0 and frame1 = ', frame0.rmsd(frame1))
        print ('rmsd frame0 and bkframe = ', frame0.rmsd(bkframe))
        print ('rmsd frame1 and bkframe = ', frame1.rmsd(bkframe))

if __name__ == "__main__":
    unittest.main()
