from pytraj import io as mdio

top = mdio.read_parm("../tests/data/Tc5b.top")

# strip all atoms but CA
top.strip_atoms("!@CA")

# save new AMBER parm to file
mdio.write_parm("./output/_Tcb5.onlyCA.top", top)
