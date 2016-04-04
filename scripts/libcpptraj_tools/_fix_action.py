from glob import glob

pxdlist = glob("Action_*.pxd")

for pxd in pxdlist:
    with open(pxd, 'r') as fh:
        lines = fh.readlines()
        for i, line in enumerate(lines):
            if line.startswith("from") and not line.startswith("from Action"):
                lines[i] = "#" + lines[i]
    with open("new/" + pxd, 'w') as fhnew:
        fhnew.writelines(lines)
