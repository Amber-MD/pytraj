
class SimplifiedResidue(object):
    def __init__(self, resname, resid, atoms, start, end):
        self.resname = resname
        self.name = self.resname
        self.resid = resid
        self.index = self.resid
        self.atoms = atoms
        self.chain = 1
        self.first_atom_idx = start
        self.last_atom_idx = end
