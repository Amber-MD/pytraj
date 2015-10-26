
class SimplifiedResidue(object):
    def __init__(self, resname, resid, atoms):
        self.resname = resname
        self.name = self.resname
        self.resid = resid
        self.index = self.resid
        self.atoms = atoms
        self.chain = 1
