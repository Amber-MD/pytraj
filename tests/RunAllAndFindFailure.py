
from subprocess import call

with open("log", 'w') as log_file, open("output.txt", 'w') as file_out:
    # get all the files starting with 'test_' and having "import unittest"
    call(['python', './get_unittest_files.py'])

    # run tests
    call(['sh','.//TestListTravis.sh'], stdout = file_out, stderr = log_file)

with open("log", 'r') as log_file:
    i_fails = 0
    for line in log_file.readlines():
        if "File" in line:
            i_fails += 1
            print (line)

print ("%s FAILs" % i_fails)
