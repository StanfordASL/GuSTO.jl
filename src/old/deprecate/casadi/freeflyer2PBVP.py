import pdb
import time
import h5py
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from casadi import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d

from setup import *
from psm import *
from auxiliary import *

def initial_guess(Xi,Xf,properties):
  R = properties['R']
  (n,m,N) = properties['n'], properties['m'], properties['N']
  (mass,J) = properties['mass'], properties['J']

  tauF = properties['tau_guess']
  tCGL,wn = clenshaw_curtis(N)
  tauCGL = t_to_tau(tCGL,0,tauF)
  
  stateguess = np.zeros((N+1,n))
  uguess = np.zeros((N+1,m))

  for (k,tau) in enumerate(tauCGL):
    stateguess[k,:] = pathPoint(tau,tauF,Xi,Xf,R)

  for (k,tau) in enumerate(tauCGL):
    if k==0 or k==N+1:
      continue
    dtau = tauCGL[k]-tauCGL[k-1]

    q = stateguess[k,6:10]
    (qx,qy,qz,qw) = q 

    qdot = (stateguess[k,6:10]-stateguess[k-1,6:10])/dtau
    Psi = 2*np.array([[qw, qz, -qy, -qx],
            [-qz, qw, qx, -qy],
            [qy, -qx, qw, -qz],
            0.5*q.T])
    stateguess[k,10:13] = np.dot(Psi,qdot)[0:3]

  for k in range(1,N+1):
    w = stateguess[k,10:13]
    dtau = tauCGL[k]-tauCGL[k-1]

    wdot = (stateguess[k,10:13]-stateguess[k-1,10:13])/dtau
    uguess[k-1,0:3] =  mass*(stateguess[k,3:6]-stateguess[k-1,3:6])/dtau
    uguess[k-1,3:6] =  np.dot(J,wdot) + np.cross(w,np.dot(J,w)) 

  return stateguess, uguess

def casadi_setup(D, wn, properties):
  # Unpack relevant properties
  N = properties['N']
  n = properties['n']
  m = properties['m']
  J = properties['J']
  (Jxx,Jyy,Jzz) = np.diag(J)
  mass = properties['mass']
  Xlb = properties['Xlb']
  Xub = properties['Xub']
  Ulb = properties['Ulb']
  Uub = properties['Uub']
  hardLimitVel = properties['hard_limit_vel']
  hardLimitAccel = properties['hard_limit_accel']
  hardLimitOmega = properties['hard_limit_omega']
  hardLimitAlpha = properties['hard_limit_omega']

  # Declare model variables
  rx = MX.sym('rx');    ry = MX.sym('ry');    rz = MX.sym('rz');
  vx = MX.sym('vx');    vy = MX.sym('vy');    vz = MX.sym('vz');
  qx = MX.sym('qx');    qy = MX.sym('qy');    qz = MX.sym('qz');    qw = MX.sym('qw');
  wx = MX.sym('wx');    wy = MX.sym('wy');    wz = MX.sym('wz');
  state = vertcat(rx,ry,rz,vx,vy,vz,qx,qy,qz,qw,wx,wy,wz)

  # tauF  = MX.sym('tauF');
  tauF = properties['tau_guess']
  tauF = 50.0 

  # Declare control variables
  Fx = MX.sym('Fx');    Fy = MX.sym('Fy');    Fz = MX.sym('Fz');
  Mx = MX.sym('Mx');    My = MX.sym('My');    Mz = MX.sym('Mz');
  control = vertcat(Fx,Fy,Fz,Mx,My,Mz)

  # System dynamics equations
  statedot = vertcat(vx,
                     vy,
                     vz,
                     1/mass*Fx,
                     1/mass*Fy,
                     1/mass*Fz,
                     0.5*(qy*wz - qz*wy + qw*wx),
                     0.5*(qz*wx - qx*wz + qw*wy),
                     0.5*(qx*wy - qy*wx + qw*wz),
                     0.5*(-qx*wx - qy*wy - qz*wz),
                     (Mx + Jyy*wy*wz - Jzz*wy*wz)/Jxx,
                     (My - Jxx*wx*wz + Jzz*wx*wz)/Jyy,
                     (Mz + Jxx*wx*wy - Jyy*wx*wy)/Jzz)

  f = Function('f', [state,control], [statedot])

  # Start with an empty NLP
  w = []    # the decision variables lbw < w < ubw
  lbw = []  # lower bound constraint on dec var
  ubw = []  # upper bound constraint on dec var
  g = []    # vector for constraints lbg < g(w,p) < ubg
  lbg = []  # lower bound on constraints
  ubg = []  # upper bound on constraints

  # Decision variables 
    # state
  X = []
  for k in range(N+1):
    Xn = MX.sym('X_'+str(k), n)
    X += [Xn]
    w += [Xn]
    lbw += Xlb.tolist()
    ubw += Xub.tolist()

    # control
  U = []
  for k in range(N+1):
    Un = MX.sym('U_' + str(k), m)
    U += [Un]
    lbw += Ulb.tolist()
    ubw += Uub.tolist()
    w   += [Un]

    # final time
  # lbw += [properties['taulb'][0]]
  # ubw += [properties['tauub'][0]]
  # w += [tauF]

  # Quaternion norm constraints
  qeps = 1e-4
  qnorm = Function('qnorm', [state], [qx**2+qy**2+qz**2+qw**2])
  for k in range(N+1):
    Xn = X[k]
    g += [qnorm(Xn)]
    lbg += [1-qeps]
    ubg += [1+qeps]

  # Speed constraints
  vel = Function('vel', [state], [vx**2+vy**2+vz**2])
  omega = Function('omega', [state], [wx**2+wy**2+wz**2])
  for k in range(N+1):
    Xn = X[k]
    g += [vel(Xn), omega(Xn)]
    lbg += [0,0]
    ubg += [hardLimitVel**2,hardLimitOmega**2]

  # Acceleration constraints 
  accel = Function('accel', [control], [Fx**2+Fy**2+Fz**2])
  alpha = Function('alpha', [control], [(Mx/Jxx)**2 + (My/Jyy)**2 + (Mz/Jzz)**2])

  for k in range(N+1):
    Un = U[k]
    g += [accel(Un), alpha(Un)]
    lbg += [0,0]
    ubg += [hardLimitAccel**2,hardLimitAlpha**2] 

  # Dynamics constraints 
  for k in range(N+1):
    # Expression for the state derivative at the collocation point
    dn = D[k,0]*X[0]
    for j in range(1,N+1):
      dn += D[k,j]*X[j]

    # Append collocation equations to constraint dynamics
    fn = f(X[k], U[k])
    g += [fn - (2.0/tauF)*dn]
    lbg += np.zeros((n)).tolist() 
    ubg += np.zeros((n)).tolist() 

  # Cost function
  energy = Function('energy', [state,control], [Fx**2+Fy**2+Fz**2+Mx**2+My**2+Mz**2])
  J = 0.0
  for k in range(N+1):
    J += wn[k]*energy(X[k],U[k])

  # Create an NLP solver
  prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
  solver = nlpsol('S', 'ipopt', prob);

  return solver, {'lbx':lbw,'ubx':ubw,'lbg':lbg,'ubg':ubg}

def casadi_solve(Xi,Xf,solver,con_bd,properties,use_guess=True):
  # Unpack problem structure
  lbx = con_bd['lbx']
  ubx = con_bd['ubx']
  lbg = con_bd['lbg']
  ubg = con_bd['ubg']

  # Unpack relevant properties
  N = properties['N']
  n = properties['n']
  m = properties['m']

  # Update BCs
  lbx[:n] = Xi.tolist();  lbx[n*N:n*(N+1)] = Xf.tolist();
  ubx[:n] = Xi.tolist();  ubx[n*N:n*(N+1)] = Xf.tolist();

  if use_guess:
    stateguess, uguess = initial_guess(Xi,Xf,properties)
    w0 = []   # initial guess of decision variables
    for k in range(N+1):
      w0 += stateguess[k].tolist()
    for k in range(N+1):
      w0 += uguess[k].tolist()

    # Solve the NLP
    soln = solver(x0=w0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
  else:
    soln = solver(lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)

  return soln, solver.stats()['success']

def trajectory(sol, properties):
  """ Parse the output from Ipopt and interpolate the polynomials """
  # Unpack relevant properties
  N = properties['N']
  n = properties['n']
  m = properties['m']

  # Parse the solution
  X_opt  = sol['x'][:n*(N+1)]
  x_opt  = X_opt[0::n]
  y_opt  = X_opt[1::n]
  z_opt  = X_opt[2::n]
  vx_opt = X_opt[3::n]
  vy_opt = X_opt[4::n]
  vz_opt = X_opt[5::n]
  qx_opt = X_opt[6::n]
  qy_opt = X_opt[7::n]
  qz_opt = X_opt[8::n]
  qw_opt = X_opt[9::n]
  wx_opt = X_opt[10::n]
  wy_opt = X_opt[11::n]
  wz_opt = X_opt[12::n]
  # tau_opt  = X_opt[-1]
  states = np.hstack((x_opt,y_opt,z_opt,vx_opt,vy_opt,vz_opt,qx_opt,qy_opt,qz_opt,qw_opt,wx_opt,wy_opt,wz_opt)) 

  U_opt = sol['x'][n*(N+1):]
  Fx = U_opt[0::m]
  Fy = U_opt[1::m]
  Fz = U_opt[2::m]
  Mx = U_opt[3::m]
  My = U_opt[4::m]
  Mz = U_opt[5::m]
  controls = np.hstack((Fx,Fy,Fz,Mx,My,Mz))

  return states, controls

def main(fn=None):
  # solve 2PBVP for graph
  if not fn:
    fn = 'V.h5' 

  with h5py.File(fn, 'r') as hf:
    V = hf['V'][:].T  # n x n_samples
  n_samples = V.shape[1]

  properties = load_params()

  # Compute N Chebyshev Gauss Lobatto points
  tCGL, wn = clenshaw_curtis(properties['N'])
  D = chebyshev_diff_matrix(properties['N'])

  # Setup solver for NLP
  solver, con_bd = casadi_setup(D, wn, properties)

  # Solve the NLP
  slew_weighting = properties['slew_weighting'] 
  max_se3_dist = properties['max_se3_dist']

  for i in range(n_samples):
    Xi = V[:,i]
    for j in range(n_samples):
      if i == j:
        continue

      Xf = V[:,j]

      # check if Xi and Xf close enough
      dist = (np.linalg.norm(Xi[0:3]-Xf[0:3]) + 
              slew_weighting * np.arccos(np.dot(Xi[6:10],Xf[6:10])))

      states = np.empty([properties['n'],0])
      controls = np.empty([properties['m'],0])
      cost = np.inf
      solved = False
      if dist < max_se3_dist:
        soln,solved = casadi_solve(Xi, Xf, solver, con_bd, properties)
        # Parse the solution
        if solved:
          cost = float(soln['f'])
          states, controls = trajectory(soln,properties)
      
      # Saving the objects
      with open('solns/soln_{}to{}.pickle'.format(i,j), 'w') as f: 
          pickle.dump([states,controls,cost,solved], f)
      print("Done with edge ({},{}){}".format(i,j,'\n'*15))

def construct_tpbvp(fn=None):
  # given nodes in a FMT solution, solve 2PBVP between each node
  if not fn:
    fn = 'V_fmt.h5' 

  with h5py.File(fn, 'r') as hf:
    V = hf['V'][:].T  # n x n_samples
  n_wps = V.shape[1]

  properties = load_params()

  # Compute N Chebyshev Gauss Lobatto points
  tCGL, wn = clenshaw_curtis(properties['N'])
  D = chebyshev_diff_matrix(properties['N'])

  # Setup solver for NLP
  solver, con_bd = casadi_setup(D, wn, properties)

  # Solve the NLP
  slew_weighting = properties['slew_weighting'] 
  max_se3_dist = properties['max_se3_dist']

  states = []
  controls = []
  run_times = []

  n,N = properties['n'],properties['N']
  eps = 0.05  # float threshold to check if BCs satisfied
  for k in range(n_wps-1):
    t0 = time.time()
    Xi,Xf = V[:,k], V[:,k+1]
    soln,solved = casadi_solve(Xi, Xf, solver, con_bd, properties)
    state, control = trajectory(soln,properties)
    states += [state]
    controls += [control]
    run_times += [time.time()-t0]


  return states, controls, run_times, solved 
