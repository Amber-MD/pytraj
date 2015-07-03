from __future__ import absolute_import
import os
from .. TrajectoryIterator import TrajectoryIterator

__all__ = ['Ala3_crd', 'Ala3_crd_top', 'tz2_ortho_nc', 'tz2_ortho_parm7']

mydir = os.path.dirname(os.path.abspath(__file__))

Ala3_crd = os.path.join(mydir, "Ala3", "Ala3.crd")
Ala3_crd_top = os.path.join(mydir, "Ala3", "Ala3.top")
tz2_ortho_nc = os.path.join(mydir, "tz2", "tz2.ortho.nc")
tz2_ortho_parm7 = os.path.join(mydir, "tz2", "tz2.ortho.parm7")
