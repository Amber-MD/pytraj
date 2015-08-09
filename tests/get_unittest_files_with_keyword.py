from glob import glob
import sys

keyword = sys.argv[1]

lines = []
testlist = glob("test_*.py")
for pyfile in testlist:
    with open(pyfile, 'r') as fh:
        txt = fh.read()
        if "import unittest" in txt and not "#import unittest" in txt:
            if keyword in txt:
                line0 = "echo ./%s \n" % pyfile
                line = "python ./%s || exit 1 \n" % pyfile
                lines.append(line0)
                lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
