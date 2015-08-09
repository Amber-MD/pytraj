from __future__ import print_function
import re
from glob import glob

pxdlist = glob("Analysis_*.pxd")
pxdlist.remove("Analysis_Matrix.pxd")
print (pxdlist)

for action_with_ext in pxdlist:
#for action_with_ext in [pxdlist[0],]:
    action = action_with_ext.split(".")[0]
    with open("Analysis_Matrix.pxd") as fmat:
        txt = fmat.read()

    txt = txt.replace("Analysis_Matrix", action)

    with open("./bk/" + action + ".pxd", 'w') as fnew:
        fnew.writelines(txt)
