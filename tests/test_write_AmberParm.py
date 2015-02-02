from pytraj import Topology
from pytraj.parms.Parm_Amber import Parm_Amber

parm = Parm_Amber()

top = Topology("data/Tc5b.top")

parm.write_parm("./output/_Tc5b.top", top)

top2 = top.copy()
top2.strip_atoms("!@CA")
parm.write_parm("./output/_Tc5b_stripbutkeepCA.top", top)
