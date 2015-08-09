__grid__ = """
TODO : add pytraj' doc
Examples (cpptraj)
-----------------
    # Grid water density around a solute.
    trajin tz2.truncoct.nc
    autoimage origin
    rms first :1-13
    # Create average of solute to view with grid.
    average avg.mol2 :1-13
    grid out.dx 20 0.5 20 0.5 20 0.5 :WAT@O

    #Generate grid from bounds command.
    trajin tz2.ortho.nc
    autoimage
    rms first :1-13&!@H= mass
    bounds :1-13 dx .5 name MyGrid out bounds.dat
    average bounds.mol2 :1-13
    # Save coordinates for second pass.
    createcrd MyCoords
    run
    # Grid using grid data set from bounds command.
    crdaction MyCoords grid bounds.xplor data MyGrid :WAT@O

    # Create non-orthogonal grid:
    trajin tz2.truncoct.nc
    reference ../tz2.truncoct.nc [REF]
    autoimage triclinic
    grid nonortho.dx boxref [REF] 50 50 50 :WAT@O pdb nonortho.pdb
"""
