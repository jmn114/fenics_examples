import numpy as np
import orthopoly as op
from pseudospectral import diffmat
from chebinterp import chebinterp as interp

def nodaltrial1d(ni,nq):
    """ Construct 1D Langrange polnomial nodal trial functions
        using the Chebyshev-Gauss-Lobatto nodes as the interpolation
        points and evaluates them at the Legendre-Gauss nodes 
 
        Input parameters:
        ni - order of polynomial interpolant space
        ng - number of Gauss points

        Outputs:
        xi - interpolation nodes
        xq - integration nodes
        wq - quadrature weights
        L  - matrix of Lagrange polynomials
        Lx - matrix of derivatives of Lagrange polynomials """

    theta = np.pi*np.arange(ni+1)/ni
    xi = -np.cos(theta) 
    a,b = op.rec_jacobi(nq,0,0)
    xq, wq = op.gauss(a,b)
    I = np.identity(ni+1)
    D = diffmat(xi)
    L = interp(xi,I,xq)
    Lx = interp(xi,D,xq)  

    return xi,xq,wq,L,Lx
