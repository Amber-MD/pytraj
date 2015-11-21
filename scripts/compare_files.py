
import os
from glob import glob

# get action files from pytraj
pyxfiles = []
with open("./pyxlist.txt", 'r') as f:
    for line in f.readlines():
        if "#" not in line and "Action_" in line:
            pyxfiles.append(line.split("\n")[0].split("/")[1])

pyxfiles = set(pyxfiles)

# get action files from cpptraj
#cppdir = os.environ['CPPTRAJHOME'] + "/src"
cppdir = os.environ['AMBERHOME'] + "/AmberTools/src/cpptraj/src/"
h_files = set([act.split("/")[-1].split(".")[0] for act in glob(cppdir + "/Action_*")])

print(h_files.difference(pyxfiles))
