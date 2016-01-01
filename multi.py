"""FEEG6002 Advance Computational Method Coursework 2015/16: PDE methods.
   Name: Alan Tan Kay Meng  Student ID: 25816322 """

import math
import numpy as np
from numpy import dot, sqrt
import scipy as sp
from scipy.linalg import block_diag
import matplotlib.pyplot as plt
np.set_printoptions(linewidth = 999999 )

# Global variable to quickly change matrix size
Mynum = 9


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
    # n = N - 2 for embed
    n = N
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
    assert np.all(CheckU == b) # working only mynum = 7 and 9 

    # Print value of the products at midpoint.
    mid = (n**2-1)/2
    print "Q1: Value of the dot product A.u1 is %5.3f at (0.5,0.5)." % (CheckU[mid])
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
print "Q2: Value of the dot product A.x2 is %5.3f at (0.5,0.5)." % (ans[(Mynum**2-1)/2])
plt.figure(2)
plt.clf()
plot_pcolor(Tfull_q2)
plt.savefig('Q2_AlanTan_25816322.pdf')
## Main program END ------------------------------------------------------------


##### Question 2 END -----------------------------------------------------------


#==============================================================================#


##### Question 3 START ---------------------------------------------------------
"""Replace the simple 4 point stencil with the following stencil in your code for
question (1) and solve PDE (1)"""

## Construct Laplacian matrix START --------------------------------------------
def get_A3(n):
    """Return matrix A for 2D Laplace equation using block diagonal
    structure, given the number of unknowns 'n' in each direction.
    """
    # Create a matrix B
    Bdiag = -60 * np.eye(n)
    Bupper1 = np.diag([16] * (n - 1), 1)
    Bupper2 = np.diag([-1] * (n - 2), 2)
    Blower1 = np.diag([16] * (n - 1), -1)
    Blower2 = np.diag([-1] * (n - 2), -2)
    B = Bdiag + Bupper1 + Blower1 + Bupper2 + Blower2

    # Creat a list [B,B,B,...,B] with n Bs
    blst = [B] * n

    # Unpack and rearrange list of Bs into diagonal of matrix A
    A = sp.linalg.block_diag(*blst)

    # Upper diagonal array offset by n: we've got (n-1) I blocks
    # each containing n ones
    Dupper1 = np.diag(16*np.ones(n * (n - 1)), n)
    Dupper2 = np.diag(-1*np.ones(n * (n - 2)), 2*n)

    # Lower diagonal array offset by -n
    Dlower1 = np.diag(16*np.ones(n * (n - 1)), -n)
    Dlower2 = np.diag(-1*np.ones(n * (n - 2)), -2*n)
    A += Dupper1 + Dlower1 + Dupper2 + Dlower2

    # Print the A matrix
    # print A 
    return A
## Construct Laplacian matrix END -----``--------------------------------------- 


## Generic solver START --------------------------------------------------------
def laplace2dq3(get_A3, get_rho, N=Mynum, Te=2):
    """Build in solver to compute value of u with determined boundary condition"""
    # Reduce the row and column of Laplacian matrix by 2 
    # Reduced row and column will be replace with embed in future
    # n = N - 2
    n = N

    # Solving for the PDE(1)
    h = 1.0/(n-1)
    A = get_A3(n) * (1/(h**2))
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
    assert np.all(CheckU == b) # work for Mynum = 7 and 9

    # Print value of the products at midpoint.
    mid = (n**2-1)/2
    print "Q3: Value of the dot product A.u3 is %5.3f at (0.5,0.5)." % (CheckU[mid])
    return Tfull
## Generic solver END ----------------------------------------------------------


## Main program START ----------------------------------------------------------
Tfull3 = laplace2dq3(get_A3, get_rho)
plt.figure(3)
plt.clf()
plot_pcolor(Tfull3)
plt.savefig('Q3_AlanTan_25816322.pdf')

## Main program END ------------------------------------------------------------


##### Question 3 END -----------------------------------------------------------


#==============================================================================#


##### Question 4 START--- ------------------------------------------------------
"""Develop a Gauss-Seidel "red-black" solver and use it to solve PDE(1)."""

def redblackA(N=Mynum):
    """Matrix A"""
    h = 1.0 / (N-1)
    Aold = get_A(N) * (1/h**2)
    Anew = np.diag(np.zeros(N**2))
    n = N**2

    for i in range(0,(n+1),2):
        for j in range(0,(n+1),2):
            Anew[i/2,j/2] = Aold[i,j]
    
    for i in range(1,(n+1),2):
        for j in range(1,(n),2):
            Anew[(i-1)/2,((n-1)/2)+((j+1)/2)] = Aold[i-1,j]

    for i in range(0,(n-1),2):
        for j in range(0,(n+1),2):
            Anew[((n+1)/2)+(i/2),j/2] = Aold[i+1,j]

    for i in range(1,(n-1),2):
        for j in range(1,(n-1),2):
            Anew[((n-1)/2)+((i+1)/2) , ((n-1)/2)+((j+1)/2)] = Aold[i,j]

    return Anew

def redblackb(N=Mynum):
    """Matrix b"""
    bold = get_rho(N)
    n = N**2
    bnew = np.zeros(n)

    for i in range(0,n,2):
        bnew[i/2] = bold[i]

    for i in range(1,n,2):
        bnew[((n-1)/2)+(i+1)/2] = bold[i]

    return bnew

def redblackb_rev(x,N=Mynum):
    n = len(x)
    bnew = np.zeros(n)

    for i in range(0,(n+1)/2):
        bnew[i*2] = x[i]

    for i in range((n+1)/2 , n):
        bnew[(i*2)-n] = x[i]

    return bnew

def iterate4(x, omega=1, N=Mynum):
    """Use the Gauss-Seidel algorithm to iterate the estimated solution
    vector x to equation A x = b, and return the improved solution.

    x : array of floats of size n
         Solution vector.
    omega : float
         Relaxation factor.

    """
    n = len(x)
    h = 1.0 / (N - 1.)
    A = redblackA(N)
    b = redblackb(N)
    
    m = (n-1)/2
    l = (n-1)
    
    for i in range(0,n):
        xsum=0
        for j in range(0,n):
            xsum =  xsum + A[i,j]*x[j] 
        xsum = xsum - A[i,i]*x[i]    
        x[i] = omega * (b[i] - xsum) / A[i,i] + (1-omega)*x[i]
   
    return x

x4, niter4, omega4 = gauss_seidel(iterate4, array, tol=1.0e-9)
x4rev = redblackb_rev(x)
Aq4 = get_A(Mynum)
h = 1. / (Mynum - 1.)
ans4 = np.dot((1/h**2)*Aq4,x4rev)

T_q4 = x4rev.reshape((Mynum, Mynum))
Tfull_q4 = embed(T_q4, 2)

for i in range(0,len(x4)):
    if ans4[i]<1e-6:
        ans4[i]=0

plt.figure(4)
plt.clf()
plot_pcolor(Tfull_q4)
plt.savefig('Q4_AlanTan_25816322.pdf')

print "Q4: Value of the dot product A.x4 is %5.3f at (0.5,0.5)." % (ans4[(Mynum**2-1)/2])
##### Question 4 END -----------------------------------------------------------