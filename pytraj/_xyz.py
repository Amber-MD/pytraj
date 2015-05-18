
class XYZ(object):
    """Wrapper of numpy array for read-only
    This class is intended to be used with Trajectory-like class

    There are limited methods compared to numpy's ndarray
    """
    def __init__(self, xyz):
        self._xyz = xyz

    def __getattr__(self, name):
        return getattr(self._xyz, name)

    def __array__(self):
        return self._xyz[:]

    def __str__(self):
        return self._xyz.__str__()

    def __repr__(self):
        return self._xyz.__str__()

    def __len__(self):
        return self._xyz.__len__()

    def __getitem__(self, idx):
        return self._xyz[idx]

    def __setitem__(self, idx, value):
        raise NotImplementedError("read only")

    def __iadd__(self, value):
        raise NotImplementedError("read only")

    def __isub__(self, value):
        raise NotImplementedError("read only")

    def __imul__(self, value):
        raise NotImplementedError("read only")
    
    def __idiv__(self, value):
        raise NotImplementedError("read only")

    def __ifloordiv__(self, value):
        raise NotImplementedError("read only")

    def __mul__(self, value):
        return self._xyz.__mul__(value)

    def __div__(self, value):
        return self._xyz.__div__(value)

    def __add__(self, value):
        return self._xyz.__add__(value)

    def __sub__(self, value):
        return self._xyz.__sub__(value)

    def __abs__(self, value):
        return self._xyz.__abs__(value)

    def __eq__(self, value):
        return self._xyz.__eq__(value)

    def __floordiv__(self, value):
        return self._xyz.__floordiv__(value)

    def __truediv__(self, value):
        return self._xyz.__floordiv__(value)
