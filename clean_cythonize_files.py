import os
# get *.pyx files
pyxfiles = []
f = open('pyxlist.txt', 'r')
try:
    for line in f.readlines():
        if "#" not in line:
            pyxfiles.append(line.split("\n")[0])
finally:
    f.close()

rootname = os.getcwd()
pytraj_home = rootname + "/pytraj/"

for ext_name in pyxfiles:
    pyxfile = pytraj_home + ext_name + ".cpp"
    #print (pyxfile)
    do_this = "mv %s trash" % pyxfile
    os.system(do_this)
