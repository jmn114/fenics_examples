from numpy import arange, outer, cos, arccos, dot
from numpy.polynomial.chebyshev import chebvander, chebval
from scipy.linalg import solve

def chebinterp(x,f,y):
    """ Interpolate a function f on a set of points x onto a set of points y
        x and y should must be in [-1,1]"""
    n = len(x)
    Tx = chebvander(x,n-1)
    a = solve(Tx,f)
    p = chebval(y,a).T
    return p
