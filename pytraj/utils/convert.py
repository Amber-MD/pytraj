
def array_to_cpptraj_range(seq):
    # use "i+1" since cpptraj use 1-based index for mask
    return ",".join((str(i+1) for i in seq))

def array_to_cpptraj_atommask(seq):
    return '@' + array_to_cpptraj_range(seq)

def array_to_cpptraj_residuemask(seq):
    return ':' + array_to_cpptraj_range(seq)

def array2d_to_cpptraj_maskgroup(arr):
    '''
    >>> array2d_to_cpptraj_maskgroup([[0, 2], [3, 4]])
    @1,3 @4,5
    '''
    import numpy as np
    arr = np.asarray(arr)
    assert arr.ndim == 2, 'must be 2d array-like'
    a0 = array_to_cpptraj_atommask(arr[0])
    a1 = array_to_cpptraj_atommask(arr[1])
    return ' '.join((a0, a1))

def atom_pairs_to_cpptraj_atommask(atom_pairs):
    """
    Notes: not validate yet

    Parameters
    ----------
    atom_pairs : 2D array-like

    Returns
    -------
    1D array of cpptraj mask

    Examples
    >>> atom_pairs_to_cpptraj_atommask([[0, 3], [5, 6]])
    ['@1 @4', '@6 @7']
    """
    import numpy as np
    if not isinstance(atom_pairs, np.ndarray):
        atom_pairs = np.asarray(atom_pairs)

    atom_pairs = (atom_pairs + 1).astype('str')

    out = np.empty(atom_pairs.shape[0], dtype=object)
    for idx, pair in enumerate(atom_pairs):
        # add "1" as 1-based indexing for cpptraj
        out[idx] = " ".join(('@'+pair[0], '@' + pair[1]))
    return out
