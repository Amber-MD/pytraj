import pytraj as pt

top = pt.load_topology("../tests/data/Tc5b.top")

# save only CA atoms
pt.write_parm("./output/_Tcb5.onlyCA.top", top["@CA"])
