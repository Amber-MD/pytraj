from glob import glob

lines = []
# turn on later
# testlist = glob("*.py") + glob("./more_compicated_examples/*.py")
testlist = glob("*.py")

# remove ./run_all_and_find_fails.py to avoid infinite loops
remove_list = ['run_examples.py', 'get_py_files.py',
               'example_calculate_chi_angle.py',
               'example_randomize_ions.py',
               'example_energy_decomposition.py',
               'example_rotate_dihedral_and_energy_calc.py',
               'example_rmsd_2trajs.py',  # do not have large files on travis
               'example_kmeans_sklearn.py',
               'example_pca_sklearn_vs_cpptraj.py',
               'example_load_file_from_url.py']

for key in remove_list:
    try:
        testlist.remove(key)
    except ValueError:
        pass

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
