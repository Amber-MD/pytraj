from pytraj.DihedralSearch import DihedralSearch
from pytraj.Topology import Topology

top = Topology("./data/Tc5b.top")

dhsearch = DihedralSearch()
dhsearch.list_known_types()
dhsearch.offset_help()

