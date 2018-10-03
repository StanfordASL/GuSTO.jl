import numpy as np
from scipy.interpolate import interp1d

def t_to_tau(t,tau0,tauf):
  """ Converts from the interval [-1,1] to the interval [tau0,tauf] """

  tau = np.zeros(len(t))
  for i in range(len(tau)):
    tau[i] = 0.5*((tauf-tau0)*t[i] + (tauf + tau0))
  return tau


def tau_to_t(tau):
  """ Converts from the interval [tau0,tauf] to the interval [-1,1] """
  t = np.zeros(len(tau))
  tau0 = tau[0]
  tauf = tau[-1]
  for i in range(len(tau)):
    t[i] = (2.0*tau[i] - (tauf-tau0))/(tauf-tau0)
  return t


def chebyshev_nodes(N):
    """ Returns N+1 CGL points on domain [1,-1] for the extrema
    of the Nth order Chebyshev polynomial. """
    x = np.zeros(N+1)
    for k in range(N+1):
        x[k] = np.cos(k*np.pi/N)
    return x


def clenshaw_curtis(N):
  """ Computes the weights for discretizing integral based on Clenshaw-Curtis
  quadrature scheme, uses CGL points """
  
  # Compute the Chebyshev-Gauss-Lobatto points
  x = chebyshev_nodes(N)
  x = x[::-1]

  # Compute optimal weights
  w = np.zeros(N+1)

  # For N even
  if N % 2 == 0:
    w[0] = 1.0/(N**2 - 1)
    w[N] = w[0]
    a = 0
  else:
    w[0] = 1.0/N**2
    w[N] = w[0]
    a = 1
  for s in range(1,(N-a)/2 + 1):
    w[s] = 2.0/N
    for j in range(1,(N-a)/2):
      w[s] += (4.0/N)*(1.0/(1.0-4.0*j**2))*np.cos(2*np.pi*j*s/N)
    w[s] += (2.0/N)*(1.0/(1-(N-a)**2))*np.cos((N-a)*s*np.pi/N)
    w[N-s] = w[s]

  return x,w


def chebyshev_diff_matrix(N): 
  """ Returns differentiation matrix D, such that D multiplied by the column
  vector [x[0],x[1],...x[N]].T gives the derivatives at each of the x points,
  and x is on the interval [-1,1] """
  x = chebyshev_nodes(N)
  D = np.zeros((N+1,N+1))
  for k in range(N+1):
    c = np.ones(N+1)
    c[0] = 2
    c[N] = 2
    for j in range(N+1):
      if k == 0 and j == 0:
        D[k,j] = (2*N**2 + 1)/6.0
      elif k == N and j == N:
        D[k,j] = -(2*N**2 + 1)/6.0
      elif k == j:
        D[k,j] = -x[k]/(2*(1-x[k]**2))
      else:
        D[k,j] = (c[k]/c[j])*(np.power(-1,j+k)/(x[k]-x[j]))

  # Modify to make the time vector [-1,1] instead of [1,-1]
  D = -D
  return D


def interpolation_matrix(N, x):
  """ Returns the matrix PHI which given a vector of node points 
  xnode = [x0 x1 x2 ..]^T, PHI*xnode will result in a vector of points
  on the interpolating polynomial at the times given in x. """
  xCGL = chebyshev_nodes(N)
  xCGL = xCGL[::-1]
  M = len(x)

  PHI = np.zeros((M,N+1))

  # Construct polynomial basis (using Lagrange polynomials as the basis)
  for i in range(N+1):
      # Construct Lagrange polynomials to get the polynomial basis at the collocation point
      phi_i = np.poly1d([1])
      for k in range(N+1):
          if k != i:
              phi_i *= (1/(xCGL[i]-xCGL[k]))*np.poly1d([1, -xCGL[k]])

      # Evalute the interpolating polynomial at the point x
      for k in range(M):
          PHI[k,i] = phi_i(x[k])

  return PHI

def poly_approximation(t,x):
  """ Computes the polynomial approximation of x(t) based on the Lagrange
  polynomial basis 
  Inputs:
  t: CGL nodes
  x: optimal points """
  N = len(t) - 1
  poly = 0
  for j in range(N+1):
    # Construct jth Lagrange polynomial
    phi_j = np.poly1d([1])
    for k in range(N+1):
      if k != j:
        phi_j *= (1/(t[j]-t[k]))*np.poly1d([1, -t[k]])
    poly += np.poly1d([x[j]])*phi_j

  return poly


def lagrange_polynomials(t):
  """ Computes the Lagrange polynomial basis functions phi(t)
  Inputs:
  t: CGL nodes """
  N = len(t) - 1
  poly = 0
  phi = []
  for j in range(N+1):
    # Construct jth Lagrange polynomial
    phi_j = np.poly1d([1])
    for k in range(N+1):
      if k != j:
        phi_j *= (1/(t[j]-t[k]))*np.poly1d([1, -t[k]])
    phi.append(phi_j)

  return phi

def polyApproximation(t,x):
  """ Computes the polynomial approximation of x(t) based on the Lagrange
  polynomial basis 
  Inputs:
  t: CGL nodes
  x: optimal points """
  N = len(t) - 1
  poly = 0
  for j in range(N+1):
    # Construct jth Lagrange polynomial
    phi_j = np.poly1d([1])
    for k in range(N+1):
      if k != j:
        phi_j *= (1/(t[j]-t[k]))*np.poly1d([1, -t[k]])
    poly += np.poly1d([x[j]])*phi_j

  return poly
