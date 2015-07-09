from pytraj.parms.ParmFile import ParmFile

# used to read Topology
# we don't use "ParmFile" class directly to avoid circular importing


class TMPParmFile(ParmFile):

    def __init__(self, *args, **kwd):
        pass
