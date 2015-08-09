from time import time


class Timer:

    def __init__(self):
        self.start = -1
        self._time_gap = -1

    def __enter__(self):
        self.start = time()
        return self

    def __exit__(self, *args):
        self._time_gap = time() - self.start

    def __call__(self, func):
        def inner(*args, **kwd):
            with self:
                func(*args, **kwd)
            print(self.time_gap())
        return inner

    def time_gap(self):
        return self._time_gap

    def __str__(self):
        return str(self._time_gap)
