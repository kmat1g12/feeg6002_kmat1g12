"""FEEG6002 Advance Computational Method Coursework 2015/16: PDE methods.
   Name: Alan Tan Kay Meng	Student ID: 25816322 """

import math
import numpy as np
import scipy as sp
from scipy.linalg import block_diag
import matplotlib.pyplot as plt

# Global variable to quickly change matrix size
Mynum = 51


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
    print "Q1: Value of the dot product A.u is %3.2f at (0.5,0.5)." % (CheckU[mid])
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
# plt.show()
# mid = (Mynum-1)/2
# print Tfull[mid,mid]
## Main program END ------------------------------------------------------------
##### Question 1 END -----------------------------------------------------------


##### Question 2 ---------------------------------------------------------------

##### Question 2 END -----------------------------------------------------------