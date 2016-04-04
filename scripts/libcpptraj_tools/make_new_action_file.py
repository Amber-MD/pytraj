# python ./make_new_action_file.py action_file_name

import sys
key = sys.argv[1]
with open("Action_Energy.pxd") as epxd, open("Action_Energy.pyx") as epyx:
    # pxd file
    tpxd = epxd.read()
    tpxd = tpxd.replace("Action_Energy", key)

    # pyx file
    tpyx = epyx.read()
    tpyx = tpyx.replace("Action_Energy", key)

    with open("./tmp/" + key + ".pxd", 'w') as newpxd, open("./tmp/" + key + ".pyx", 'w') as newpyx:
        newpxd.write((tpxd))
        newpyx.write((tpyx))
