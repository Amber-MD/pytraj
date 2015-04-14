from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing

has_parmed = has_("chemistry")

if not has_parmed:
    PytrajWarningMissing("`chemistry`")

def load_ParmEd(parm_name):
    if has_parmed:
        from chemistry import load_file
        return load_file(parm_name)
    else:
        return None
