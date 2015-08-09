from pytraj.base import *


class TestFrame(Frame):
    pass


frame = TestFrame()

setattr(frame, 'top', Topology())

frame.top

frame.top = Topology("./data/Tc5b.top")
frame.top.summary()
