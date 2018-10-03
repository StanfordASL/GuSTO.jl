import numpy as np

def quat_multiply(q1,q2):
  # Eq. 169-171 in Shuster
  # q1 * q2 <==> R(q1) * R(q2)

  q_product = np.zeros((4))
  q_product[0:3] = q1[3]*q2[0:3] + q2[3]*q1[0:3] - np.cross(q1[0:3],q2[0:3])
  q_product[3] = q1[3]*q2[3] - np.dot(q1[0:3],q2[0:3])
  return q_product

def quat_exp(v):
  # gives quaternion q corresponding to orientation obtained
  # from an initial orientation [0 0 0 1] by rotating of an angle norm(v)
  # around fixed direction v/norm(v)
  # Eq. 9 in "Time-Optimal Reorientation of a Spacecraft Using an Inverse Dynamics Optimization Method"
  phi = np.linalg.norm(v)
  sin_half_phi = np.sin(0.5*phi)
  if phi == 0:
    return np.array([0.,0,0,1])
  return np.append(sin_half_phi*v/phi,np.cos(0.5*phi))

def quat_inv(q):
  # Eq. 177 in Shuster
  qinv = np.copy(q)
  qinv[0:3]*=-1
  return qinv/np.linalg.norm(qinv)

def quat2angle(q):
  return 2*np.arctan2(np.linalg.norm(q[0:3]), q[3])

def quat2vec(q,canonicalize=False):
  # canonicalize: true corresponds to positive scalar component q for dual representation
  if np.linalg.norm(q[0:3]) == 0:
    return zeros((3))

  a = np.sign(q[4])*q[0:3]/np.linalg.norm(q[0:3]) if canonicalize else q[0:3]/np.linalg.norm(q[0:3])
  return a

def quat_log(q):
  # gives vector having direction of the Euler's axis and magnitude equal
  # to Euler's angle of the orientation from [0 0 0 1] to q
  # Eq. 9 in "Time-Optimal Reorientation of a Spacecraft Using an Inverse Dynamics Optimization Method"

  if np.linalg.norm(q[0:3]) == 0:
    return np.zeros((3))
  return quat2angle(q)*quat2vec(q)

def quat_interp(q0,q1,t):
  # p371 in "A General Construction Scheme for Unit Quaternion Curves with Simple High Order Derivatives"
  # NOTE: q1 * q2 <==> R(q1) * R(q2)

  if (t>1 or t<0):
    print("Supplied time t must lie in [0,1]!\n")
    return q0 if t<0 else q1
  
  quat_error = quat_multiply(quat_inv(q0),q1)
  q_out = quat_multiply(q0, quat_exp(t*quat_log(quat_error)))
  return q_out
