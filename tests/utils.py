import os

def fn(name):
    # return absolute dir of ./data/name
    return os.path.dirname(__file__) + '/data/' + name
