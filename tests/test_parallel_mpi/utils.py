import os


def fn(name):
    # Note: this function use different path compared to
    # the one in ../utils.py
    # return absolute dir of ./data/name
    return os.path.dirname(__file__) + '/../data/' + name


tc5b_trajin = fn('Tc5b.x')
tc5b_top = fn('Tc5b.top')
tz2_trajin = fn('tz2.nc')
tz2_top = fn('tz2.parm7')
tz2_ortho_trajin = fn('tz2.ortho.nc')
tz2_ortho_top = fn('tz2.ortho.parm7')
