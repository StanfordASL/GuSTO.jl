import numpy as np
from quat_functions import *

def pathPoint(t,tau,xi,xf,R):
  # t: time to query
  # tau: final time
  x = np.copy(xi) 
  x[0] = xf[0] + xf[3]*(t - tau) + (((4*R[0,0]*(xf[3] - xi[3]))/tau + (6*R[0,0]*(xi[0] - xf[0] + tau*xi[3]))/tau**2)*(t - tau)**2)/(2*R[0,0]) + (((6*R[0,0]*(xf[3] - xi[3]))/tau**2 + (12*R[0,0]*(xi[0] - xf[0] + tau*xi[3]))/tau**3)*(t - tau)**3)/(6*R[0,0])
  x[1] = xf[1] + xf[4]*(t - tau) + (((4*R[1,1]*(xf[4] - xi[4]))/tau + (6*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau**2)*(t - tau)**2)/(2*R[1,1]) + (((6*R[1,1]*(xf[4] - xi[4]))/tau**2 + (12*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau**3)*(t - tau)**3)/(6*R[1,1])
  x[2] = xf[2] + xf[5]*(t - tau) + (((4*R[2,2]*(xf[5] - xi[5]))/tau + (6*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau**2)*(t - tau)**2)/(2*R[2,2]) + (((6*R[2,2]*(xf[5] - xi[5]))/tau**2 + (12*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau**3)*(t - tau)**3)/(6*R[2,2])
  x[3] = xf[3] + (((4*R[0,0]*(xf[3] - xi[3]))/tau + (6*R[0,0]*(xi[0] - xf[0] + tau*xi[3]))/tau**2)*(t - tau))/R[0,0] + (((6*R[0,0]*(xf[3] - xi[3]))/tau**2 + (12*R[0,0]*(xi[0] - xf[0] + tau*xi[3]))/tau**3)*(t - tau)**2)/(2*R[0,0])
  x[4] = xf[4] + (((4*R[1,1]*(xf[4] - xi[4]))/tau + (6*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau**2)*(t - tau))/R[1,1] + (((6*R[1,1]*(xf[4] - xi[4]))/tau**2 + (12*R[1,1]*(xi[1] - xf[1] + tau*xi[4]))/tau**3)*(t - tau)**2)/(2*R[1,1])
  x[5] = xf[5] + (((4*R[2,2]*(xf[5] - xi[5]))/tau + (6*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau**2)*(t - tau))/R[2,2] + (((6*R[2,2]*(xf[5] - xi[5]))/tau**2 + (12*R[2,2]*(xi[2] - xf[2] + tau*xi[5]))/tau**3)*(t - tau)**2)/(2*R[2,2])

  # TODO(acauligi): Check to make sure slerp function is correct and matches Eigen library/Ames convention
  x[6:10] = quat_interp(xi[6:10],xf[6:10],t/tau)
  return x
