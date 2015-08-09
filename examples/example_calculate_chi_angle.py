import pytraj as pt

# load DNA structure
# calculate glycosidic torsion angle
trajin = "../tests/data/Test_NAstruct/adh026.3.crd"
parm = "../tests/data/Test_NAstruct/adh026.3.top"

traj = pt.iterload(trajin, parm)

print(pt.calc_chin(traj))
print(pt.calc_chin(traj, resrange=range(3, 21, 3)))
