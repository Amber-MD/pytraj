from glob import glob

for pyfile in glob("*.py"):
    with open(pyfile, 'r') as fh:
        print("python ./%s" % pyfile)
