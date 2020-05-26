import pytraj as pt
from memory_profiler import profile

orig_top = pt.load_topology('../../tests/data/tz2.parm7')


def test():
    for i in range(100000000):
        print(i)
        top = orig_top.copy()
        top.strip(":1")


test()
