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

    >>> # for data in pjob.compute(): print(data)
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
