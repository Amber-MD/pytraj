from __future__ import absolute_import 
from pytraj import AtomSelect
import numpy as np
from .gdt import gdt

def calc_score(frame0=None, frame1=None, mask="*", 
               top=None, score="gdtscore", *args, **kwd):
    #print "API will be changed"
    if score == 'gdtscore':
        _score = 1
    elif score == 'tmscore':
        _score = 2
    elif score == 'maxsub':
        _score = 3

    ats = AtomSelect(top=top)
    ats.selected_frame = frame0
    arr0 = np.asarray(ats.select(mask).flatten(), dtype=np.float32)

    ats.selected_frame = frame1
    arr1 = np.asarray(ats.select(mask).flatten(), dtype=np.float32)
    return gdt(arr0, arr1, 1, len(arr0)/3, _score)[0]
