from glob import glob


lines = []
testlist = glob("test_*.py")
for pyfile in testlist:
    with open(pyfile, 'r') as fh:
        txt = fh.read()
        if "import unittest" in txt and not "#import unittest" in txt:
            line0 = "echo ./%s \n" % pyfile
            line = "python ./%s \n" % pyfile
            lines.append(line0)
            lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
