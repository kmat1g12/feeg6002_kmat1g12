"""FEEG6002 Advance Computational Method Coursework 2015/16: PDE methods.
   Name: Alan Tan Kay Meng	Student ID: 25816322 """

import math
import numpy as np
import scipy as sp
from scipy.linalg import block_diag
import matplotlib.pyplot as plt


# Question 1
""" Solve PDE(1) using 2D Laplace equation solver."""

# Generic solver

def embed(T, Te=2):
    """Embed the array T giving the temperature at the inner nodes in
    the domain into a larger array including the boundary temperatures
    """
    N = T.shape[0] + 2
    Tfull = np.zeros((N,N))
    # Tfull[((N/2)-1)][((N/2)-1)] = Te
    Tfull[1:-1, 1:-1] = T
    return Tfull

def laplace2d(get_A, get_rho, N=50, Te=2):
    """Solve the Laplace equation on a 2D grid, with T=0 at all
    boundaries except y=0, where T=Te, and return an 2D array of size
    NxN giving the temperature distribution throughout the domain.
    """
    n = N - 2
    A = get_A(n)
    b = get_rho(n, Te)
    U = sp.linalg.solve(A, b)
    T = U.reshape((n, n))
    # Tfull = embed(T, Te)
    # return Tfull
    return T

# Plotting functions

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

# Setup the A matrix for the system of linear equations

def get_A(n):
    """Return matrix A for 2D Laplace equation using block diagonal
    structure, given the number of unknowns 'n' in each direction.
    """
    Bdiag = -4 * np.eye(n)
    Bupper = np.diag([1] * (n - 1), 1)
    Blower = np.diag([1] * (n - 1), -1)
    B = Bdiag + Bupper + Blower
    # Creat a list [B,B,B,...,B] with n Bs
    blst = [B] * n
    # Unpack the list of diagonal blocks 'blst'
    # since block_diag expects each block to be passed as separate
    # arguments. This is the same as doing block_diag(B,B,B,...,B)
    A = sp.linalg.block_diag(*blst)
    # Upper diagonal array offset by n: we've got (n-1) I blocks
    # each containing n ones
    Dupper = np.diag(np.ones(n * (n - 1)), n)
    # Lower diagonal array offset by -n
    Dlower = np.diag(np.ones(n * (n - 1)), -n)
    A += Dupper + Dlower
    return A

# Setup the b vector for the system of linear equations

def get_rho(n, Te=2):
    """Return column vector of size n^2 containing the boundary conditions."""
    b = np.zeros(n**2)
    b[:n] = -Te
    return b

# Solve the equation and plot

Tfull = laplace2d(get_A, get_rho)
plt.figure(1)
plt.clf()
plot_pcolor(Tfull)
#plt.savefig('fig06-01.pdf')
plt.show()
print Tfull
