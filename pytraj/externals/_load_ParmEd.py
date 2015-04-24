from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing

def load_ParmEd(parm_name):
    has_parmed = has_("chemistry")
    if has_parmed:
        from chemistry import load_file
        return load_file(parm_name)
    else:
        if not has_parmed:
            PytrajWarningMissing("`chemistry`")
        return None
