import pytraj as pt

traj = pt.iterload_remd(
    "../tests/data/Test_RemdTraj/rem.nc.000",
    top="../tests/data/Test_RemdTraj/ala2.99sb.mbondi2.parm7",
    T="300.0")

# cpptraj will automatically find frames with 300K
print(traj.temperatures)
