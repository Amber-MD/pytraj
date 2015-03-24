class PytrajError(Exception):
    def __init__(self, msg="Pytraj Error"):
        self.msg = msg

    def __str__(self):
        return self.msg

class PathError(PytrajError):
    pass

class EmptyTopologyError(PytrajError):
    pass
