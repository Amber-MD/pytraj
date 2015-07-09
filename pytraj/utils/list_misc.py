def flatten_list(arr0):
    """return flattened list. Mostly for testing purpose"""
    # modified from
    # http://goo.gl/etJsvx (stackoverflow)
    import itertools
    return list(itertools.chain.from_iterable(arr0))
