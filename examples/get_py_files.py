from glob import glob


lines = []
testlist = glob("*.py")
# remove ./run_all_and_find_fails.py to avoid infinite loops
remove_list = ['run_all_and_find_fails.py', 'get_py_files.py', 
               'dihedral.py',]

for key in remove_list:
    testlist.remove(key)

# turn off mpi test
for fname in testlist:
    if "mpi" in fname:
        testlist.remove(fname)

for pyfile in testlist:
    line = "python ./%s \n" % pyfile
    lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
