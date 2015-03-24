# distutils: language = c++

"shared methods betwween FrameArray and Trajin classes"""
from pytraj.trajs.Trajout import Trajout
from pytraj.six_2 import set

def _savetraj(self, filename="", fmt='unknown', overwrite=False):
    if fmt == 'unknown':
        # convert to "UNKNOWN_TRAJ"
        fmt = fmt.upper() + "_TRAJ"
    else:
        fmt = fmt.upper()

    with Trajout(filename=filename, top=self.top, fmt=fmt, 
                 overwrite=overwrite, more_args=None) as trajout:
        for idx, frame in enumerate(self):
            trajout.writeframe(idx, frame, self.top)

def _get_temperature_set(self):
    return set(self.temperatures) 
