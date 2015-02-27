from glob import glob


lines = []
testlist = glob("*.py")
# remove ./run_all_and_find_fails.py to avoid infinite loops
testlist.remove("run_all_and_find_fails.py")
testlist.remove("get_py_files.py")

for pyfile in testlist:
    line = "python ./%s \n" % pyfile
    lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
