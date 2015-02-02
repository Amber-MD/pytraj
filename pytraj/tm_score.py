# adapted from CSB package
# dir: http://csb.codeplex.com/
# turn-off numpy.sum, numpy.power
# use builtin `sum`
# TODO : add more acknowledgement
from pytraj.DistRoutines import distance_frames

def _tm_d0(Lmin):
    if Lmin > 15:
        d0 = 1.24 * (Lmin - 15.0)**(1.0 / 3.0) - 1.8
    else:
        d0 = 0.5
    return max(0.5, d0)

def tm_score(x, y, L=None, 
             d0=None, mask=None, 
             top=None,
             mask1=None, mask2=None,
             top1=None, top2=None):
    # TODO : create new top if `top` is string
    """
    Evaluate the TM-score of two Frames (no fitting is done).
    Return : float
    
    Parameters:
    ----------
    x : Frame instance
        Frame 1
    y : Frame instance 
        Frame 2
    L: length for normalization (default: ?)
    d0: d0 in Angstroms (default: calculate from C{L})
    """
    import numpy as np

    if mask is not None:
        mask1 = mask
        mask2 = mask
        top1 = top
        top2 = top
        xx = x.strip_atoms(mask1, top1, copy=True)
        yy = y.strip_atoms(mask2, top2, copy=True)
    else:
        have_top1_2 = (top1 != None) & (top2 != None)
        have_mask1_2 = (mask1 != None) & (mask2 != None)
        if have_top1_2 and have_mask1_2:
            xx = x.strip_atoms(mask1, top1, copy=True)
            yy = y.strip_atoms(mask2, top2, copy=True)
        else:
            xx = x
            yy = y
    if not L:
        L = x.n_atoms
    if not d0:
        d0 = _tm_d0(L)
    # convert Python array to numpy array 
    d = np.array(distance_frames(xx, yy))
    return sum(1 / (1 + (d / d0) ** 2)) / L
