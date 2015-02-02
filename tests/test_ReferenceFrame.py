import os
from pytraj.AtomMask import AtomMask
from pytraj.ReferenceFrame import ReferenceFrame
from pytraj.Topology import Topology
from pytraj.FileName import FileName
from pytraj.Frame import Frame
from pytraj.CpptrajState import CpptrajState
from pytraj.Energy import Energy_Amber
from pytraj.ArgList import ArgList

topname = "./data/Tc5b.top"
filename = "./data/md1_prod.Tc5b.x"

ref = ReferenceFrame()
