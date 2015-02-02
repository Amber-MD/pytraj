class A (object):
    def __init__(self, x):
        self._x = x

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

if __name__ == "__main__":
    a = A(3)
    assert a.x == a._x
    a.x = 4
    assert a.x == a._x
    a._x = 5
    assert a.x == a._x
