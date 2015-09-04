from .Frame cimport _Frame, Frame

cpdef list pairwise_rmsd(traj, bint fit=True):
    '''cpptraj has this method but I re-write to avoid memory overhead
    (due to I don't know yet how to use TrajectoryIterator with cpptraj without 
    copying)
    '''
    cdef Frame frame0, frame1
    cdef int n_frames = traj.n_frames
    cdef list data

    data = []

    for idx, frame0 in enumerate(traj(stop=n_frames-1)):
        for frame1 in traj.iterframe(start=idx+1):
            if fit:
                data.append(frame0.rmsd(frame1))
            else:
                data.append(frame0.rmsd_nofit(frame1))
    return data
