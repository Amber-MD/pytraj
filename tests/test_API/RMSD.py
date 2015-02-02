import numpy as np
from pytraj.Topology import Topology
from pytraj.AtomMask import AtomMask
from pytraj.ArgList import ArgList
from pytraj.Frame import Frame
from pytraj.Trajin_Single import Trajin_Single
from pytraj.ReferenceFrame import ReferenceFrame

def _cpptraj_rmsd(topname, refname, trajname, first_nframes=100, mask='', use_mass=True, fit=True):
    """
    DEMO API
    TODO: update this doc
    Output: array of rmsd values
    
    Input:
    =====
    topname : Amber topology file name
    refname : Reference file name
    trajname : trajectory file name
    first_nframes : calc rmsd for first_nframes
    mask : AtomMask string
    use_mass : bool
    fit : fit or no_fit

    """
    # not yet supported atom mask
    # atm = AtomMask(mask)
    # make topology instance
    top = Topology(topname)

    # create reference frame
    ref = ReferenceFrame()
    ref.load_ref(refname, top)
    refframe = ref.frame

    # load trajectory file
    trajin = Trajin_Single()
    trajin.load(trajname, ArgList(), top)
    trajin.prepare_for_read(True)

    #create frame to store data from traj iteration
    frame = Frame()
    frame.set_frame_v(top, trajin.has_vel(), trajin.n_repdim)

    # do rmsd
    nframes = first_nframes if first_nframes < trajin.total_frames else trajin.total_frames
    arr = np.empty(nframes)

    trajin.begin_traj(False)
    trajin.print_info_line()
    for i in range(nframes):
        # not added "mask" yet
        # just for demo
        trajin.get_next_frame(frame)
        if fit:
            arr[i] = frame.rmsd(refframe, use_mass)
        else:
            arr[i] = frame.rmsd_no_fit(refframe, use_mass)
    trajin.end_traj()
    return arr

def _mdtraj_rmsd(topname, refname, trajname, first_nframes=100, mask='', use_mass=False, fit=True):
    # TODO: how to load *.crd file?
    import mdtraj as md
    pass
    #top = md.load_prmtop(topname)
    #ref = '' 
    #traj = md.load(trajname, top)
    #return md.rmsd(traj, ref, nframe=first_nframes)

def rmsd(topname, refname, trajname, first_nframes=100, mask='', 
         use_mass=False, fit=True, program='cpptraj'):
    if program == 'cpptraj':
        return _cpptraj_rmsd(topname, refname, trajname, first_nframes, use_mass, fit)
    elif program == 'mdtraj':
        return _mdtraj_rmsd(topname, refname, trajname, first_nframes, use_mass, fit)


if __name__ == '__main__':
    topname = "../data/Tc5b.top"
    refname = "../data/Tc5b.nat.crd"
    trajname = "../data/md1_prod.Tc5b.x"
    #mask = ":3-18@CA"
    n_frames = 10
    arr = rmsd(topname, refname, trajname, first_nframes=n_frames, use_mass=True, program='cpptraj')
    # load first 100 rmsd values from cpptraj calculation
    rmsd_cpptraj = np.loadtxt("../data/rmsd_allatoms.dat").transpose()[1][:n_frames]
    # STATUS: can not reproduce cpptraj's calculation
    # %error ~ 3%
    print("cpptraj: ", rmsd_cpptraj)
    print("pytraj: ", arr)
    print("difference: ", arr - rmsd_cpptraj)
