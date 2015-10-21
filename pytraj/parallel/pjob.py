class PJob(object):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.parallel import PJob
    >>> traj = pt.load_sample_data('tz2')
    >>> tasklist = []
    >>> tasklist.append((pt.radgyr, traj))
    >>> tasklist.append((pt.molsurf, traj, '@CA'))

    >>> # perform each action on each CPUs (total 2 CPUs)
    >>> pjob = PJob(tasklist)

    >>> for data in pjob.compute(): print(data)
    (0, array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
            18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
            18.91669565,  18.87069722]))
    (1, array([ 458.51409637,  459.64784573,  456.54690793,  467.72939574,
            462.45908781,  458.70327554,  454.40514806,  455.15015576,
            468.70566447,  456.0058624 ]))
    '''
    def __init__(self, tasklist):
        self.tasklist = tasklist

    def _worker(self, rank):
        my_task = self.tasklist[rank]
        func = my_task[0]
        return (rank, func(*my_task[1:]))

    def compute(self):
        from multiprocessing import Pool

        n_cores = len(self.tasklist)
        p = Pool(len(self.tasklist))

        result = p.map(self._worker, [rank for rank in range(n_cores)])
        return result
