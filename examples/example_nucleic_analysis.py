import pytraj as pt
fn = "../tests/data/Test_NAstruct/adh026.3.pdb"
traj = pt.load(fn, fn)

d = pt.nucleic_acid_analysis(traj, dtype='dataset')
print (d)

# get major groove
print (d.grep("major"))

# get major groove mean
print (d.grep("major").mean())

# get bot major and minir grooves
print (d.grep(['major', 'minor']).mean())
