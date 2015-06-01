from glob import glob


lines = []
testlist = glob("test_*.py")
base_line = "import unittest"

for pyfile in testlist:
    with open(pyfile, 'r') as fh:
        txt = fh.read()
        # only test if file having keyword `import unittest`
        # but exclude ones with `#import unittest`
        # or one with "import unittest # pragma no_test"
        if base_line in txt:
            if not base_line + " # pragma no_test" in txt and not "#" + base_line in txt:
                print (pyfile)
                line0 = "echo ./%s \n" % pyfile
                line = "python ./%s \n" % pyfile
                lines.append(line0)
                lines.append(line)

with open("./TestListTravis.sh", 'w') as fh:
    fh.writelines(lines)
