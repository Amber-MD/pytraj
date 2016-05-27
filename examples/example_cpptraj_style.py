import pytraj as pt

refname  = '../tests/data/Tc5b.nat.inpcrd'
fn = '../tests/data/Tc5b.x'
prmtop  = '../tests/data/Tc5b.prmtop'

traj = pt.iterload(fn, top=prmtop)
ref = pt.iterload(refname, top=prmtop)

cm = """
rms reference @CA
distance :3 :7
angle :3 :8 :12
multidihedral phi psi
molsurf @CA
"""

# need to specify ref=ref
# TODO: reflist
data = pt.compute(cm, traj, ref=ref)

print(data)
