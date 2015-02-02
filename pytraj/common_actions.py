"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd, do_fitting
>>> from pytraj.common_actions import translate
"""
from pytraj.DistRoutines import distance 
from pytraj.gdt.calc_score import calc_score

calc_distance = distance
calc_protein_score = calc_score

## calc
# calc_distance
# calc_radgyr(frame, mask, top)
# calc_dihedral(a1,a2,a3,a4)
# clalc_angle

## do
# do_fitting(frame, ref=Frame())
# do_translation
# do_rotation
# 

