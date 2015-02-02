"""This is a thin wrapper of Trajin_Single
We need to sub-class Trajin_Single to use FrameArray
(we called Trajin_Single from FrameArray, so we can not call FrameArray back from 
Trajin_Single)
"""
from pytraj.Trajin_Single import Trajin_Single
#from pytraj.FrameArray import FrameArray

class TrajReadOnly(Trajin_Single):
    def __init__(self, *args, **kwd):
        pass
