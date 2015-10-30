# require pytraj, pysander and parmed

import pytraj as pt
import numpy as np

traj = pt.iterload("../tests/data/Ala3/Ala3.crd",
                   "../tests/data/Ala3/Ala3.top")

print(traj.n_atoms, traj.top.n_residues)
print(pt.multidihedral(traj).to_dict())

t0 = traj[:1]
deg_ene = []

flist = []

try:
    import sander
    import parmed
    # scan from -180 to 180, every 5 deg
    # calculate dihedral energy for each conformation
    #
    for deg in range(-180, 180, 5):
        pt.rotate_dihedral(t0, "custom:3:omega:" + str(deg))

        flist.append(t0[0].copy())

        en = pt.energy_decomposition(t0,
                                     igb=8,
                                     verbose=False,
                                     parm=traj.top.filename)['dihedral'][0]
        print(deg, en)
        deg_ene.append((deg, en))

    arr = np.array(deg_ene).T

    print(len(flist))

    # write multiple pdb to a single file to visualize
    # rotation with vmd.
    pt.write_traj("output/test.pdb", flist,
                  top=traj.top,
                  overwrite=True,
                  options='model')

    # plot the PMF
    try:
        import matplotlib
        from pytraj import plot
        arr[1] = arr[1] - np.min(arr[1])
        pt.plot(0, 1, data=arr)
        pt.savefig("output/test.png")
    except ImportError:
        pass
except ImportError:
    pass
