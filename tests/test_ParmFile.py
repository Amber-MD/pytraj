import unittest
from pytraj.parms.ParmFile import ParmFile
from pytraj.Topology import Topology
from pytraj.decorators import not_yet_supported
from pytraj import io as mdio


def test_generator(top):
    for atom in top.atom_iter():
        pass
        # print atom

    print("\nres_iterator for top")
    for res in top.residue_iter():
        pass
        # print res

    print("\nmol_iterator for top")
    for mol in top.mol_iter():
        pass
        # print mol


top = Topology()
parm = ParmFile()
filename = "./data/Tc5b.top"
parm.readparm(filename, top)
assert top.n_atoms == 304
parm.help()


class TestParmFile(unittest.TestCase):
    def test_readpdb(self):
        print("test_readpdb +++++++++++++++++++++++++++++++")
        traj0 = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        mdio.write_traj(filename="./output/test_0_after.pdb",
                        traj=traj0[0],
                        top=traj0.top,
                        overwrite=True)
        pdbtop = mdio.read_parm("./output/test_0_after.pdb")
        assert pdbtop.n_atoms == 304
        print(pdbtop.n_atoms)
        top2 = pdbtop.copy()
        top2.strip_atoms("!@CA")
        print(top2.n_atoms)
        assert top2.n_atoms == 20

        # make sure that making Topology instance from pdb file is FINE
        traj0 = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        traj1 = mdio.iterload(
            "./data/md1_prod.Tc5b.x", "./output/test_0_after.pdb")
        assert traj0[0].coords == traj1[0].coords
        assert traj0[0].coords == traj1[0].coords

    def test_0(self):
        top = Topology()
        parm = ParmFile()
        print(parm.formats())

        filename = "./data/Tc5b.top"

        print("Read topology file and dump to 'top'")
        parm.readparm(filename, top)
        assert top.n_atoms == 304

        print("\natom_iterator for top")

        print(parm.filename())

    def test_write_CHARMMPSF(self):
        print("write parm for CHARMM")
        parm.writeparm(top=top, filename="output/test.psf", format="CHARMMPSF")
        top2 = Topology()
        parm.readparm(filename="output/test.psf", top=top2)
        assert top2.n_atoms == 304
        test_generator(top2)

    def test_write_AMBERPARM(self):
        top2 = Topology()
        parm.writeparm(top=top, filename="output/test.parm7")
        parm.readparm(filename="output/test.parm7", top=top2)
        assert top2.n_atoms == 304
        test_generator(top2)

    #@not_yet_supported
    def test_write_GROMACS(self):
        top2 = Topology()
        parm.writeparm(top=top, filename="output/test.gmxtop")
        #parm.readparm(filename="output/test.gmxtop", top=top2)
        #assert top2.n_atoms == 304
        # test_generator(top2)

        #@not_yet_supported
    def test_write_mol2(self):
        top2 = Topology()
        parm.writeparm(top=top, filename="output/test.mol2")
        #parm.readparm(filename="output/test.mol2", top=top2)
        #assert top2.n_atoms == 304
        # test_generator(top2)

        #@not_yet_supported
    def test_write_TINKER(self):
        top2 = Topology()
        parm.writeparm(top=top, filename="output/test.arc")
        #parm.readparm(filename="output/test.arc", top=top2)
        #assert top2.n_atoms == 304
        # test_generator(top2)


if __name__ == "__main__":
    unittest.main()
