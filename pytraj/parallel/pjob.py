class PJob(object):
    def __init__(self, tasklist):
        self.tasklist = tasklist

    def _worker(self, rank):
        my_task = self.tasklist[rank]
        func = my_task[0]
        return (rank, func(*my_task[1:]))

    def compute(self):
        '''
        '''
        from multiprocessing import Pool

        n_cores = len(self.tasklist)
        p = Pool(len(self.tasklist))

        result = p.map(self._worker, [rank for rank in range(n_cores)])
        return result
