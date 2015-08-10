Misc
====

**read Gaussian output to `pytraj.Trajectory`**::
    
    import pytraj as pt
    pt.tools.read_gaussian_output('./gau.log')
    pt.tools.read_gaussian_output('./gau.log', top='test.pdb')
