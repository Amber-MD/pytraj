
def array_to_cpptraj_range(seq):
    return ",".join((str(i+1) for i in seq))

def array_to_cpptraj_atommask(seq):
    return '@' + array_to_cpptraj_range(seq)
