Misc
====

**read Gaussian output to `pytraj.Trajectory`**::
    
    import pytraj as pt
    pt.tools.read_gaussian_output('./gau.log')
    pt.tools.read_gaussian_output('./gau.log', top='test.pdb')

**plotting**::

    from pytraj import plot
    plot(0, 1, data=data) # same as plt.plot(data[0], data[1])
