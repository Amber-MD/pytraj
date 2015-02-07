from time import time

class Timer:
    def __init__(self):
        self.start = -1
        self.time_gap = -1

    def __enter__(self):
        self.start = time()
        return self

    def __exit__(self, *args):
        self.time_gap = time() - self.start
