from glob import glob


lines = []
testlist = glob("test_*.py")
base_line = "import unittest"

for pyfile in testlist:
    with open(pyfile, 'r') as fh:
        txt = fh.read()
        if not base_line + " # pragma no_test" in txt and not "#" + base_line in txt:
            line0 = "echo ./%s \n" % pyfile
            line = "python ./%s \n" % pyfile
            lines.append(line0)
            lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
