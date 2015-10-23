import pytraj as pt

try:
    import mdtraj as md

    # load mdtraj object
    m_traj = md.load_netcdf('../tests/data/tz2.ortho.nc',
                            '../tests/data/tz2.ortho.parm7')
    print(m_traj)

    # convert to pytraj object
    # you can use a pdb file, a mol2 file, ... as Topology too
    # as long as pytraj/cpptraj supports
    traj = pt.Trajectory(xyz=m_traj.xyz.astype('f8'), top='../tests/data/tz2.ortho.parm7')
    print(traj)

    # perform 'action' on traj
    traj.autoimage()

    # copy data back to mdtraj object
    m_traj.xyz = traj.xyz[:]

    # mdtraj has very fast rmsd calculation, you can pass pytraj'traj object
    # to mdtraj to 'borrow' its action too.
    # note that pytraj/cpptraj use Angstrom for unit while mdtraj use nm
    print(md.rmsd(traj, traj, 0))

except ImportError:
    print("does not have mdtraj")
