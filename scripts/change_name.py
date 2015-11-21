# change *_pxd.pxd to *.pxd

import os
from glob import glob
import re

for fname in glob("*_pxd.pxd"):
    fnewname = re.sub("_pxd", "", fname)
    # print fnewname
    print("change %s to %s" % (fname, fnewname))
    os.rename(fname, fnewname)
