import numpy as np

def load_params():
  # problem parameters
  n = 6 
  m = 3 

  N = 100

  R = np.eye(3) # cost matrix for 2PBVP

  slew_weighting = 1.2 
  max_se2_dist = 1.50 

  # robot parameters
  # $ASTROBEE/freeflyer/astrobee/config/worlds/granite.config
  J = np.array([[0.1083, 0.0, 0.0],
                [0.0, 0.1083, 0.0],
                [0.0, 0.0, 0.1083]])
  Jinv = np.linalg.inv(J)
  Jzz = J[2,2]
  mass = 14.4 

  hard_limit_vel   = 0.50
  hard_limit_accel = 0.10
  hard_limit_omega = 45*np.pi/180
  hard_limit_alpha = 50*np.pi/180

  # state box constraints
  Xlb = np.array([-np.inf,-np.inf,-np.inf,-hard_limit_vel/np.sqrt(2),-hard_limit_vel/np.sqrt(2),-hard_limit_omega])
  Xub = -1*Xlb

  # control box constraints
  Ulb = np.array([-mass*hard_limit_accel/np.sqrt(2), -mass*hard_limit_accel/np.sqrt(2), -Jzz*hard_limit_alpha])
  Uub = -1*Ulb 

  # time constraints
  taulb = np.array([0.1])
  tauub = np.array([100.0])
  tau_guess = 600.0 

  return {'n':n, 'm':m, 'N':N, 'R':R, 'slew_weighting':slew_weighting, 'max_se2_dist':max_se2_dist, 'J':J, 'Jinv':Jinv, 'mass':mass, 
          'hard_limit_vel':hard_limit_vel, 'hard_limit_accel':hard_limit_accel, 'hard_limit_omega':hard_limit_omega, 'hard_limit_alpha':hard_limit_alpha, 
          'Xlb':Xlb, 'Xub':Xub, 'Ulb':Ulb, 'Uub':Uub, 'taulb':taulb, 'tauub':tauub, 'tau_guess':tau_guess}
