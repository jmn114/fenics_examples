import numpy as np
from pylab import plot, figure, show
from nodaltrial1d import *
from scipy.linalg import cho_factor, cho_solve


if __name__ == '__main__':

    # Polynomial order of Lagrange trial/test functions   
    p = 50

    # Number of time steps
    N = 1000

    # Final time
    T = 1

    # Time step
    dt = float(T)/N

    # Uniform time steps
    t = np.linspace(0,T,N+1)

    # IMEX Runge-Kutta parameters
    mu  = (2-np.sqrt(2))/2
    lam = 1 - 1/(2*mu)

    # Diffusion coefficient
    dcoeff = 0.1

    # Interpolation points, quadrature rule, Lagrange trial functions 
    # and their derivatives in one dimension
    xi,xq,wq,L,Lx = nodaltrial1d(p,p+10)

    # Mass matrix
    M = np.dot(L.T*wq,L)

    # Stiffness matrix
    K = dcoeff*np.dot(Lx.T*wq,Lx)

    # Matrix to be inverted at every time step and intermediate stage
    A = M + mu*dt*K

    # Cholesky factor for matrix
    R = cho_factor(A)

    # Left "inflow"
    am = lambda t: t*(T-t)

    # Right "outflow"
    ap = lambda t: t*(T-t)

    # Nonlinearity
    f = lambda u: -u**2/2

    # Use constant initial data
    u = np.ones((p+1,N+1))

    fig = figure(1)
    ax = fig.add_subplot(1,1,1)

    # Do time stepping
    for n in range(N):

        # Evaluate nonlinearity at old time
        fni = f(u[:,n])
        fn = np.dot(Lx.T*wq,np.dot(L,fni))

        Mun = np.dot(M,u[:,n])

        # right hand side
        rhs = Mun + mu*dt*fn

        that = t[n]+mu*dt

        # Correct boundary terms
        rhs[0]  += mu*dt*( fni[0]  - am(that) )
        rhs[-1] -= mu*dt*( fni[-1] - ap(that) )

        # Compute intermediate value of solution
        uhat = cho_solve(R,rhs)

        # Update nonlinearity
        fhati = f(uhat)
        fhat = np.dot(Lx.T*wq,np.dot(L,fhati))

        # Make a new right hand side  
        rhs = Mun - dt*(1-mu)*np.dot(K,uhat) + dt*(lam*(fn-fhat)+fhat)

        # Correct boundary terms
        rhs[0]  += dt*( (1-lam)*fhati[0]  + lam*fni[0]  - \
                        mu*am(t[n+1]) - (1-mu)*am(that) )
        rhs[-1] -= dt*( (1-lam)*fhati[-1] + lam*fni[-1] - \
                        mu*ap(t[n+1]) - (1-mu)*ap(that) )

        u[:,n+1] = cho_solve(R,rhs)


    ax.plot(xi,u)
    show()