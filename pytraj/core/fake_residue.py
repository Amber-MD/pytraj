
class SimplifiedResidue(object):
    '''SimplifiedResidue, only has info about name, index, start and end atom index

    Examples
    --------
    >>> import pytraj as pt
    >>> top = pt.datafiles.load_tz2_ortho().top
    >>> for atom in top.atoms:
    ...     print(atom.name, atom.residue.name, atom.residue.index)

    '''
    def __init__(self, resname, resid, atoms, start, end):
        self.resname = resname
        self.name = self.resname
        self.resid = resid
        self.index = self.resid
        self.atoms = atoms
        self.chain = 1
        self.first_atom_idx = start
        self.last_atom_idx = end
