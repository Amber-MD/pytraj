from pytraj.utils import has_
from pytraj.utils import _import
from pytraj.warnings import PytrajWarningMissing

_, chem = _import("chemistry")

def indices(self):
    # https://github.com/ParmEd/ParmEd/issues/117#issuecomment-96871159
    index_array = []
    for term in self:
        index_array.append([])
        for attr in dir(term):
            if attr.startswith('atom'):
                index_array[-1].append(getattr(term, attr).idx)
    return index_array # You can cast to numpy array first if you want

# Turn it into a descriptor on the TrackedList type
try:
    chem.TrackedList.indices = property(indices)
    del indices
except:
    pass

def load_ParmEd(parm_name):
    has_parmed = has_("chemistry")
    if has_parmed:
        return chem.load_file(parm_name)
    else:
        if not has_parmed:
            PytrajWarningMissing("`chemistry`")
        return None
