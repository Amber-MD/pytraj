import pytraj as pt

# turn on this to see cpptraj's stdout
#pt.set_cpptraj_verbose(True)

# use iterload for memory saving
traj = pt.iterload("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")

# write multiple pdbs to a single file
# just like multiple pdbs when you load from rcsb
# you can use this pdb file to view in VMD without loading
# prmtop/gro/psf/...

pt.write_traj("./output/test.pdb", traj,
              format='pdb',
              mode='model',
              overwrite=True)
