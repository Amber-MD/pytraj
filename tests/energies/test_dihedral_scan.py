import pytraj as pt
import numpy as np

traj = pt.iterload("./data/Ala3/Ala3.crd",
                   "./data/Ala3/Ala3.top")

print (traj.n_atoms, traj.top.n_residues)
print (pt.multidihedral(traj).to_dict())

t0 = traj[:1]
deg_ene = []

flist = []

for deg in range(-180, 180, 5):
    pt.rotate_dihedral(t0, "custom:3:omega:" + str(deg))

    flist.append(t0[0].copy())

    en = pt.energy_decomposition(t0, igb=8, verbose=False, parm=traj.top.filename)['dihedral'][0]
    print (deg, en)
    deg_ene.append((deg, en))

arr = np.array(deg_ene).T

print (len(flist))
pt.write_traj("test.pdb", flist, top=traj.top, 
              overwrite=True,
              mode='model')

arr[1] = arr[1] - np.min(arr[1])
pt.plotting.plot(arr[0], arr[1])
pt.plotting.plt.show()
