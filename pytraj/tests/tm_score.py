# from: https://searchcode.com/codesearch/raw/28101488/
# turn-off numpy.sum, numpy.power
# use builtin `sum`

def _tm_d0(Lmin):
    
    #from numpy import power
     
    if Lmin > 15:
        #d0 = 1.24 * power(Lmin - 15.0, 1.0 / 3.0) - 1.8
        d0 = 1.24 * (Lmin - 15.0)**(1.0 / 3.0) - 1.8
    else:
        d0 = 0.5
    return max(0.5, d0)

def tm_score(x, y, L=None, d0=None):
    """
    Evaluate the TM-score of two conformations as they are (no fitting is done).
    
    @param x: 3 x n input array
    @type x: numpy array
    @param y: 3 x n input array
    @type y: numpy array
    @param L: length for normalization (default: C{len(x)})
    @type L: int
    @param d0: d0 in Angstroms (default: calculate from C{L})
    @type d0: float
    
    @return: computed TM-score
    @rtype: float
    """
    #from numpy import sum

    if not L:
        L = len(x)
    if not d0:
        d0 = _tm_d0(L)
    d = distance(x, y)

    return sum(1 / (1 + (d / d0) ** 2)) / L
