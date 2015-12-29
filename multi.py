"""FEEG6002 Advance Computational Method Coursework 2015/16: PDE methods.
   Name: Alan Tan Kay Meng	Student ID: 25816322 """

import math
import numpy as np
from numpy import dot, sqrt
import scipy as sp
from scipy.linalg import block_diag
import matplotlib.pyplot as plt

# Global variable to quickly change matrix size
Mynum = 25


##### Question 1 ---------------------------------------------------------------
""" Solve PDE(1) using 2D Laplace equation solver."""

## Construct Laplacian matrix START --------------------------------------------
def get_A(n):
    """Return matrix A for 2D Laplace equation using block diagonal
    structure, given the number of unknowns 'n' in each direction.
    """
    # Create a matrix B
    Bdiag = -4 * np.eye(n)
    Bupper = np.diag([1] * (n - 1), 1)
    Blower = np.diag([1] * (n - 1), -1)
    B = Bdiag + Bupper + Blower

    # Creat a list [B,B,B,...,B] with n Bs
    blst = [B] * n

    # Unpack and rearrange list of Bs into diagonal of matrix A
    A = sp.linalg.block_diag(*blst)

    # Upper diagonal array offset by n: we've got (n-1) I blocks
    # each containing n ones
    Dupper = np.diag(np.ones(n * (n - 1)), n)

    # Lower diagonal array offset by -n
    Dlower = np.diag(np.ones(n * (n - 1)), -n)
    A += Dupper + Dlower
    return A

def get_rho(n, Te=2):
    """Return column vector of size n^2 containing the boundary conditions."""
    # Create a list of zeros for vector b
    b = np.zeros(n**2)

    # Apply boundary condition, midpoint of list b is 2
    b[(((n**2)-1))/2] = Te
    return b
## Construct Laplacian matrix END ----------------------------------------------


## Generic solver START --------------------------------------------------------
def laplace2d(get_A, get_rho, N=Mynum, Te=2):
    """Build in solver to compute value of u with determined boundary condition"""
    # Reduce the row and column of Laplacian matrix by 2 
    # Reduced row and column will be replace with embed in future
    n = N - 2

    # Solving for the PDE(1)
    h = 1.0/(n-1)
    A = get_A(n) * (1/(h**2))
    b = get_rho(n, Te)
    U = sp.linalg.solve(A, b)

    # Reshape the u vector into nxn matrix for heat map plotting
    T = U.reshape((n, n))

    # Embed the surrounding of U matrix into zeros
    Tfull = embed(T, Te)

    # Verify that dot function of A matrix and U vector
    # return the same rho value at midpoint
    CheckU = np.dot(A,U)

    # Filter very small value into zeros
    for i in range(0,len(CheckU)):
        if (abs(CheckU[i]) < 1e-12):
            CheckU[i] = 0

    # Validate that product of A and U matrix is the same as rho vector
    # Will give warning if it is not the same
    assert np.all(CheckU == b)

    # Print value of the products at midpoint.
    mid = (n**2-1)/2
    print "Q1: Value of the dot product A.u is %5.3f at (0.5,0.5)." % (CheckU[mid])
    return Tfull

def embed(T, Te=2):
    """Embed the surrounding of the Laplacian matrix with zeros"""
    N = T.shape[0] + 2

    # Create a matrix with NxN of zeros 
    # Append the target matrix into the center of the zeros matrix
    Tfull = np.zeros((N,N))
    Tfull[1:-1, 1:-1] = T
    return Tfull
## Generic solver END ----------------------------------------------------------


## Plotting functions START ----------------------------------------------------
def plot_pcolor(Tfull):
    """Plot temperature in the domain using pcolor"""
    N = Tfull.shape[0]
    x = y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x,y)
    plt.pcolor(X, Y, Tfull)
    plt.axis('scaled')
    plt.colorbar()
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('T(x,y) on %dx%d grid' % (N,N))

def plot_wireframe(Tfull):
    """Plot temperature in the domain using plot_wireframe"""
    from mpl_toolkits.mplot3d import axes3d
    N = Tfull.shape[0]
    x = y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x,y)
    # Construct and return a function suitable for interactive demo
    def plot(elev=25, azim=50):
        fig = plt.figure(1, figsize=(14, 8))
        plt.clf()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Tfull)
        ax.view_init(elev=elev, azim=azim)
        plt.axis('scaled')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.title('T(x,y) on %dx%d grid' % (N,N))
    plot()
    return plot
## Plotting functions END ------------------------------------------------------


## Main program START ----------------------------------------------------------
# Solve the equation and plot
# Code execution happens here
Tfull = laplace2d(get_A, get_rho)
plt.figure(1)
plt.clf()
plot_pcolor(Tfull)
plt.savefig('Q1_AlanTan_25816322.pdf')
## Main program END ------------------------------------------------------------
##### Question 1 END -----------------------------------------------------------

#==============================================================================#

##### Question 2 ---------------------------------------------------------------
"""Instead of using a built-in solver method, implement and use a successive 
over-relaxation solver"""

## Iterative Solver START ------------------------------------------------------

def iterate(x, omega=1, N=Mynum):
    """Use the Gauss-Seidel algorithm to iterate the estimated solution
    vector x to equation A x = b, and return the improved solution.

    x : array of floats of size n
         Solution vector.
    omega : float
         Relaxation factor.

    """
    n = len(x)
    h = 1.0 / (N - 1.)
    A = (1/h**2)*get_A(N)
    
    m = (n-1)/2
    l = (n-1)
    
    x[0] = omega * -( A[0,1]*x[1] + A[0,N]*x[N] ) / A[0,0] + (1-omega)*x[0]

    for i in range(1,N):
        x[i] = omega * -( A[i,i-1]*x[i-1] + A[i,i+1]*x[i+1] + A[i,i+N]*x[i+N] ) / A[i,i] + (1-omega)*x[i]

    for i in range(N, m):
        x[i] = omega * -( A[i,i-N]*x[i-N] + A[i,i-1]*x[i-1] + A[i,i+1]*x[i+1] + A[i,i+N]*x[i+N] ) / A[i,i] + (1-omega)*x[i]

    x[m] = omega * ( 2 -( A[m,m-N]*x[m-N] + A[m,m-1]*x[m-1] + A[m,m+1]*x[m+1] + A[m,m+N]*x[m+N] ) ) / A[m,m] + (1-omega)*x[m]

    for i in range(m+1, n-N):
        x[i] = omega * -( A[i,i-N]*x[i-N] + A[i,i-1]*x[i-1] + A[i,i+1]*x[i+1] + A[i,i+N]*x[i+N] ) / A[i,i] + (1-omega)*x[i]

    for i in range(n-N,l):
        x[i] = omega * -( A[i,i-1]*x[i-1] + A[i,i+1]*x[i+1] + A[i,i-N]*x[i-N] ) / A[i,i] + (1-omega)*x[i]

    x[l] = omega * -( A[l,l-1]*x[l-1] + A[l,l-N]*x[l-N] ) / A[l,l] + (1-omega)*x[l]

    return x

def gauss_seidel(iterate, x, tol=1.0e-9, relaxation=True):
    """ x, niter, omega = gauss_seidel(iterate, x, tol=1.0e-9, omega=1.0)

    Gauss-Seidel method for solving [A]{x} = {b}.

    The matrix [A] should be sparse. User must supply the
    function iterate(x, omega) that returns the improved {x},
    given the current {x}. 'omega' is the relaxation factor.
    """
    omega = 1.
    k = 10
    p = 1
    for i in range(1,501):
        xold = x.copy()
        x = iterate(x, omega)
        dx = sqrt(dot(x - xold, x - xold))
        if dx < tol:
            return x, i, omega
        if relaxation:
            # Compute of relaxation factor after k+p iterations
            if i == k:
                dx1 = dx
            if i == k + p:
                dx2 = dx
                omega = 2.0 / (1.0 + sqrt(1.0 - (dx2 / dx1)**(1.0 / p)))
    print 'Gauss-Seidel failed to converge'

## Iterative Solver END --------------------------------------------------------

## Main program START ----------------------------------------------------------
array = np.ones(Mynum**2)
x, niter, omega = gauss_seidel(iterate, array, tol=1.0e-9)
Aq2 = get_A(Mynum)
h = 1. / (Mynum - 1.)
ans = np.dot((1/h**2)*Aq2,x)

T_q2 = x.reshape((Mynum, Mynum))
Tfull_q2 = embed(T_q2, 2)

for i in range(0,len(ans)):
    if ans[i]<1e-6:
        ans[i]=0

print "Q2: Value of the dot product A.x is %5.3f at (0.5,0.5)." % (ans[(Mynum**2-1)/2])
plt.figure(2)
plt.clf()
plot_pcolor(Tfull_q2)
plt.savefig('Q2_AlanTan_25816322.pdf')
## Main program END ------------------------------------------------------------
##### Question 2 END -----------------------------------------------------------

#==============================================================================#

##### Question 3 START ---------------------------------------------------------

##### Question 3 END -----------------------------------------------------------