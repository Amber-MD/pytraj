"""having common actions such as rmsd, fitting, ...
>>> from pytraj.common_actions import calc_rmsd, do_fitting
>>> from pytraj.common_actions import translate
"""
from functools import partial
from pytraj import adict
from pytraj.misc import calculate, calc_dssp
from pytraj.DistRoutines import distance 
from pytraj.gdt.calc_score import calc_score

calc_distance = partial(calculate, 'distance')
calc_dih = partial(calculate, 'dihedral')
calc_dihedral = calc_dih
calc_radgyr = partial(calculate, 'radgyr')
calc_angle = partial(calculate, 'angle')
calc_molsurf = partial(calculate, 'molsurf')
calc_molsurf = partial(calculate, 'molsurf')
calc_distrmsd = partial(calculate, 'distrmsd')
calc_watershell = partial(calculate, 'watershell')
calc_protein_score = calc_score
