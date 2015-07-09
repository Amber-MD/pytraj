
def array_to_cpptraj_range(seq):
    # use "i+1" since cpptraj use 1-based index for mask
    return ",".join((str(i+1) for i in seq))

def array_to_cpptraj_atommask(seq):
    return '@' + array_to_cpptraj_range(seq)
