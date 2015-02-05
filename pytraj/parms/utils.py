from pytraj.parms.ParmFile import ParmFile
from pytraj.Topology import Topology
from pytraj.cpptraj_dict import ParmFormatDict
from pytraj.externals.six import string_types

parmfile = ParmFile()
parmtypes = ParmFormatDict.keys()

def readparm(filename):
    """return Topology instance"""
    top = Topology()
    parmfile.readparm(top, filename)
    return top

def writeparm(filename="", parm=None, fmt="AMBER", *args):
    """
    Parameters:
    ----------
    filename : str, output filename, default=""
    parm : str or Topology instance, default=None
    fmt : parm format, default="AMBER"
    *args : optional args (not supported yet)
    """
    if isinstance(parm, string_types):
        top = Topology(parm)
    elif isinstance(parm, Topology):
        top = parm

    # use uppercase 
    fmt = fmt.upper()
    # make sure having `fmt`
    if not parmtypes.has_key(fmt):
        # trying to add "FILE" to `fmt` 
        fmt += "FILE"

    # convert to "enum"
    fmtid  = ParmFormatDict[fmt]
    parmfile.writeparm(top, filename, fmtid, 0)

