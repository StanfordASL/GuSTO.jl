import numpy as np

def load_params():
  # problem parameters
  n = 13
  m = 6

  N = 50

  R = np.eye(3) # cost matrix for 2PBVP

  slew_weighting = 1.2 
  max_se3_dist = 1.50 

  # robot parameters
  J = np.array([[0.1083, 0.0, 0.0],
                [0.0, 0.1083, 0.0],
                [0.0, 0.0, 0.1083]])
  Jinv = np.linalg.inv(J)
  mass = 7.0

  hard_limit_vel   = 0.50
  hard_limit_accel = 0.10
  hard_limit_omega = 45*np.pi/180
  hard_limit_alpha = 50*np.pi/180

  # state box constraints
  Xlb = np.array([-np.inf,-np.inf,-np.inf, 
  -hard_limit_vel/np.sqrt(3),-hard_limit_vel/np.sqrt(3),-hard_limit_vel/np.sqrt(3),
  -1.0,-1.0,-1.0,0.0,
  -hard_limit_omega/np.sqrt(3),-hard_limit_omega/np.sqrt(3),-hard_limit_omega/np.sqrt(3)])
  Xub = np.array([np.inf,np.inf,np.inf, 
  hard_limit_vel/np.sqrt(3),hard_limit_vel/np.sqrt(3),hard_limit_vel/np.sqrt(3),
  1.0,1.0,1.0,1.0,
  hard_limit_omega/np.sqrt(3),hard_limit_omega/np.sqrt(3),hard_limit_omega/np.sqrt(3)])

  # control box constraints
  Jmin = np.min(np.diag(J))
  Ulb = np.array([-mass*hard_limit_accel/np.sqrt(3), -mass*hard_limit_accel/np.sqrt(3), -mass*hard_limit_accel/np.sqrt(3), -Jmin*hard_limit_alpha/np.sqrt(3), -Jmin*hard_limit_alpha/np.sqrt(3), -Jmin*hard_limit_alpha/np.sqrt(3)])
  Uub = np.array([mass*hard_limit_accel/np.sqrt(3), mass*hard_limit_accel/np.sqrt(3), mass*hard_limit_accel/np.sqrt(3),  Jmin*hard_limit_alpha/np.sqrt(3), Jmin*hard_limit_alpha/np.sqrt(3), Jmin*hard_limit_alpha/np.sqrt(3)])

  # time constraints
  taulb = np.array([0.1])
  tauub = np.array([100.0])
  tau_guess = 50.0 

  return {'n':n, 'm':m, 'N':N, 'R':R, 'slew_weighting':slew_weighting, 'max_se3_dist':max_se3_dist, 'J':J, 'Jinv':Jinv, 'mass':mass, 
          'hard_limit_vel':hard_limit_vel, 'hard_limit_accel':hard_limit_accel, 'hard_limit_omega':hard_limit_omega, 'hard_limit_alpha':hard_limit_alpha, 
          'Xlb':Xlb, 'Xub':Xub, 'Ulb':Ulb, 'Uub':Uub, 'taulb':taulb, 'tauub':tauub, 'tau_guess':tau_guess}
