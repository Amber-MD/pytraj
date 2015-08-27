from glob import glob

lines = []
# turn on later
#testlist = glob("*.py") + glob("./more_compicated_examples/*.py")
testlist = glob("*.py")

# remove ./run_all_and_find_fails.py to avoid infinite loops
remove_list = ['run_all_and_find_fails.py', 'get_py_files.py',
               'example_calculate_chi_angle.py',
               'example_randomize_ions.py',
               'example_load_file_from_url.py']

for key in remove_list:
    testlist.remove(key)

# turn off mpi test
for fname in testlist:
    if "mpi" in fname:
        testlist.remove(fname)
    else:
        with open(fname, 'r') as fh:
            if "# tag: no travis test" in fh.readline():
                testlist.remove(fname)

for pyfile in testlist:
    line = "python ./%s \n" % pyfile
    lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
