### CONVENTIONS ###
# quaternion: q = [rho q4] where rho is unit vector
# MRP: stereographic projection for MRP uses domain S^3 \ {(0,0,0,-1)}

# References
# A Survey of Attitude Representations; Shuster 1993
# Spacecraft Dynamics and Control; de Ruiter, Damaren, & Forbes 2013
# Time-Optimal Reorientation of a Spacecraft Using an Inverse Dynamics Optimization Method; Boyarko, Romano and Yakimenko 2011
# A General Construction Scheme for Unit Quaternion Curves with Simple High Order Derivatives; Kim et al. 1995
# Modified Rodrigues Parameters: An Efficient Representation of Orientation in 3D VIsion and Graphics; Terzakis et al. 2018

function skew(v::Vector{T}) where T
  # generates 3x3 skew symmetric matrix for vector
  return [0     -v[3]   v[2];
          v[3]   0      -v[1];
         -v[2]   v[1]   0]
end

function dcm2angle(dcm::Matrix{T}) where T
  # 1.3.2.2 on p21 in Spacecraft Dynamics and Control by Anton de Ruiter
  # Eq. 102 in Shuster
  return acos(0.5*(trace(dcm) - 1))
end


function dcm2vec(dcm::Matrix{T}) where T
  # 1.3.2.2 on p21 in Spacecraft Dynamics and Control by Anton de Ruiter
  # Eq. 103a in Shuster
  phi = dcm2angle(dcm)

  if abs(phi^2) < 0.00001 
    # vector undefined if phi == 0, but physically meaningless
    return zeros(T,3)
  elseif abs(phi) == pi
    # three possibilities if phi == pi, pick one
    # Eq. 104 in Shuster
    # Eq. 1.29 in SDAC_deRuiter on p21

    a1 = sqrt(0.5*(dcm[1,1] + 1))
    a2 = sqrt(0.5*(dcm[2,2] + 1))
    a3 = sqrt(0.5*(dcm[3,3] + 1))

    # a1 > 0
    a = [abs(a1); sign(dcm[1,2])*abs(a2); sign(dcm[1,3])*abs(a3)]
    # a2 > 0
    [sign(dcm[1,2])*abs(a1); abs(a2); sign(dcm[2,3])*abs(a3)]
    # a3 > 0
    [sign(dcm[1,3])*abs(a1); sign(dcm[2,3])*abs(a2); abs(a3)]
    return a./norm(a)
  end

  a = zeros(T,3)
  a[1] = 0.5*(dcm[2,3] - dcm[3,2])/sin(phi)
  a[2] = 0.5*(dcm[3,1] - dcm[1,3])/sin(phi)
  a[3] = 0.5*(dcm[1,2] - dcm[2,1])/sin(phi)
  return a
end


function dcm2quat(dcm::Matrix{T}, canonicalize::Bool=false) where T
  # canonicalize: true corresponds to positive scalar component q for dual representation
  # Eq. 163-168 in Shuster

  phi,axis = dcm2angle(dcm),dcm2vec(dcm)

  q = [sin(0.5*phi).*axis; cos(0.5*phi)]
  q = canonicalize ? sign(q[4])*q./norm(q) : q./norm(q)
  return q
end


function dcm2mrp(dcm::Matrix{T}, canonicalize::Bool=false) where T
  return quat2mrp(dcm2quat(dcm, canonicalize))
end


function vec2quat(axis::Vector{T}, phi::T, canonicalize::Bool=false) where T
  q = [sin(0.5*phi)*axis/norm(axis); cos(0.5*phi)]
  q = canonicalize ? sign(q[4])*q./norm(q) : q./norm(q)
  return q
end


function quat2angle(q::Vector{T}) where T
  return 2*atan2(norm(q[1:3]), q[4])
end


function quat2vec(q::Vector{T}, canonicalize::Bool=false) where T
  # canonicalize: true corresponds to positive scalar component q for dual representation
  if norm(q[1:3])^2 < 0.00001
    return zeros(T,3)
  end

  a = canonicalize ? sign(q[4])*q[1:3]/norm(q[1:3]) : q[1:3]/norm(q[1:3])
  return a
end


function quat2dcm(q::Vector{T}, attitude_convention::String="orientation") where T
  # Eq. 97 in Shuster (note sin term sign flipped based on skew matrix convention)
  # Eq. 1.26 on p20 in Spacecraft Dynamics and Control by Anton de Ruiter
  phi = quat2angle(q)
  n = quat2vec(q)
  dcm = cos(phi)*eye(3) + (1-cos(phi))*n*n' - sin(phi)*skew(n)

  if attitude_convention == "rotation"
    dcm = dcm'
  elseif attitude_convention != "orientation"
    @warn "Must specify if attitude is rotation or orientation! Defaulting to orientation."
  end
  return dcm
end


function vec2dcm(e::Vector{T}, phi::T) where T
  q = vec2quat(e, phi)
  return quat2dcm(q)
end


function quat_derivative(q::Vector{T}, w::Vector{T}) where T
  # Eq. 305 in Shuster
  omega_skew = [-skew(w) w; -w' 0]
  return 0.5*omega_skew*q
end


function quat_multiply(q1::Vector{T}, q2::Vector{T}) where T
  # Eq. 169-171 in Shuster
  # q1 * q2 <==> R(q1) * R(q2)

  q_product = zeros(T,4)
  q_product[1:3] = q1[4]*q2[1:3] + q2[4]*q1[1:3] - cross(Vector(q1[1:3]), Vector(q2[1:3]))
  q_product[4] = q1[4]*q2[4] - dot(q1[1:3], q2[1:3])
  return q_product
end


function quat_inv(q::Vector{T}) where T
  # Eq. 177 in Shuster
  qinv = copy(q)
  qinv[1:3]*=-1
  return qinv
end


function quat_error(q1::Vector{T}, q2::Vector{T}) where T
  # corresponds to  q2 = quat_multiply(dq,q1)
  # quat_multiply(q2, quat_inv(q1))
  
  # Spacecraft Dynamics and Control - Sidi
  #prod_matrix = [q1[4] q1[3] -q1[2] q1[1]; -q1[3] q1[4] q1[1] q1[2]; q1[2] -q1[1] q1[4] q1[3]; -q1[1] -q1[2] -q1[3] q1[4]]
  #dq = prod_matrix * [-q2[1:3]; q2[4]]

  # Space Vehicle Dynamics and Control - Wie
  prod_matrix = [q1[4] q1[3] -q1[2] -q1[1]; -q1[3] q1[4] q1[1] -q1[2]; q1[2] -q1[1] q1[4] -q1[3]; q1[1] q1[2] q1[3] q1[4]]
  dq = prod_matrix * q2
  return dq
end


function quat_rotate(q::Vector{T}, v::Vector{T}, attitude_convention::String="orientation") where T
  # Eq. 183 in Shuster
  # v: 3 element vector to be rotated

  qr = q
  if attitude_convention == "rotation"
    qr = quat_inv(q)
  elseif attitude_convention != "orientation"
    @warn "Must specify if quaternion represents rotation or orientation! Defaulting to orientation."
  end

  v_bar = [v; 0]
  return quat_multiply(qr, quat_multiply(v_bar, quat_inv(qr)))[1:3]
end


function quat_interp(q0::Vector{T}, q1::Vector{T}, t::T) where T
  # p371 in "A General Construction Scheme for Unit Quaternion Curves with Simple High Order Derivatives"
  # NOTE: q1 * q2 <==> R(q1) * R(q2)

  if (t>1 || t<0)
    print("Supplied time t must lie in [0,1]!\n")
    return t<0 ? q0 : q1
  end
  
  qerr = quat_multiply(quat_inv(q0),q1)
  return quat_multiply(q0, quat_exp(t*quat_log(qerr)))
end


function quat_exp(v::Vector{T}) where T
  # gives quaternion q corresponding to orientation obtained
  # from an initial orientation [0 0 0 1] by rotating of an angle norm(v)
  # around fixed direction v/norm(v)
  # Eq. 9 in "Time-Optimal Reorientation of a Spacecraft Using an Inverse Dynamics Optimization Method"
  phi = norm(v)
  if phi^2 < 0.00001
    return Vector{T}([0;0;0;1])
  end
  return Vector{T}([sin(0.5*phi)*v/phi;cos(0.5*phi)])
end


function quat_log(q::Vector{T}) where T
  # gives vector having direction of the Euler's axis and magnitude equal
  # to Euler's angle of the orientation from [0 0 0 1] to q
  # Eq. 9 in "Time-Optimal Reorientation of a Spacecraft Using an Inverse Dynamics Optimization Method"

  if norm(q[1:3])^2 < 0.00001
    return zeros(T,3)
  end
  return quat2angle(q)*quat2vec(q)
end


function quat2mrp(q::Vector{T}) where T
  # Eq. 249 in Shuster
  # Eq. 22 in Terzakis et al.
  if abs(q[4]+1) < 0.000001
    @warn "Stereographic projection undefined for (0,0,0,-1)"
    return -q[1:3]/(1-q[4])
  end
  return q[1:3]/(1+q[4])
end


function mrp2quat(p::Vector{T}) where T
  # Eq. 25 & 26 in Terzakis et al.
  s = (1-norm(p)^2)/(1+norm(p)^2)
  v = 2*p/(1+norm(p)^2)
  return [v;s]
end


function mrp2dcm(p::Vector{T}, attitude_convention::String="orientation") where T
  # Eq. 255b in Shuster
  # NOTE (6/26): switched +4 term to -4 term because that's what I
  # found everywhere else and unit tests pass with this only

  dcm = eye(3) + (8*skew(p)^2 - 4*(1 - norm(p)^2)*skew(p))/(1+norm(p)^2)^2

  if attitude_convention == "rotation"
    dcm = dcm'
  elseif attitude_convention != "orientation"
    @warn "Must specify if attitude is rotation or orientation! Defaulting to orientation."
  end
  return dcm
end


function mrp_derivative(p::Vector{T}, w::Vector{T}) where T
  # Eq. 338 in Shuster (in terms of the body-referenced angular velocity)
  return 0.25*((1-norm(p)^2)*w - 2*cross(w,p) + 2*dot(w,p)*p)
end


function mrp_multiply(p1::Vector{T}, p2::Vector{T}) where T
  # Eq. 257 in Shuster
  # p1 * p2 <==> R(p1) * R(p2)
  # TODO(acauligi): check if denominator is close to 0. for edge case

  return ((1-norm(p2)^2)*p1 + (1-norm(p1)^2)*p2 - 2*cross(p1,p2)) /
    (1+norm(p1)^2*norm(p2)^2 - 2*dot(p1,p2))
end

# function Aq(J::Matrix{T},X::Vector{T})
#   Jxx,Jyy,Jzz = diag(J)
#   qx,qy,qz,qw = X[1:4]
#   wx,wy,wz = X[5:7]
#   Jq = 0.5*[0 wz -wy wx qw -qz qy;
#           -wz 0 wx wy qz qw -qx;
#            wy -wx 0 wz -qy qx qw;
#           -wx -wy -wz 0 -qx -qy -qz]
#   Jw = [0 0 0 0 0                     (Jyy-Jzz)*wz/Jxx   (Jyy-Jzz)*wy/Jxx;
#         0 0 0 0 -(Jxx-Jzz)*wz/Jyy  0                     -(Jxx-Jzz)*wx/Jyy;
#         0 0 0 0 (Jxx-Jyy)*wy/Jzz   (Jxx-Jyy)*wx/Jzz                     0]
#   return [Jq;Jw]  
# end

# function Apb(J::Matrix{T},X::Vector{T})
#   Jxx,Jyy,Jzz = diag(J)
#   px,py,pz = X[1:3]
#   wx,wy,wz = X[4:6]
#   Jpb = [(px*wx)/2+(py*wy)/2+(pz*wz)/2 wz/2+(px*wy)/2-(py*wx)/2 (px*wz)/2-wy/2-(pz*wx)/2 px^2/4-py^2/4-pz^2/4+1/4 (px*py)/2-pz/2 py/2+(px*pz)/2;
#         (py*wx)/2-(px*wy)/2-wz/2 (px*wx)/2+(py*wy)/2+(pz*wz)/2 wx/2+(py*wz)/2-(pz*wy)/2 pz/2+(px*py)/2 -px^2/4+py^2/4-pz^2/4+1/4 (py*pz)/2-px/2;
#         wy/2-(px*wz)/2+(pz*wx)/2 (pz*wy)/2-(py*wz)/2-wx/2 (px*wx)/2+(py*wy)/2+(pz*wz)/2 (px*pz)/2-py/2 px/2+(py*pz)/2 -px^2/4-py^2/4+pz^2/4+1/4]
#   Jw = [0 0 0 0 (Jyy*wz-Jzz*wz)/Jxx (Jyy*wy-Jzz*wy)/Jxx;
#         0 0 0 -(Jxx*wz-Jzz*wz)/Jyy 0 -(Jxx*wx-Jzz*wx)/Jyy;
#         0 0 0 (Jxx*wy-Jyy*wy)/Jzz (Jxx*wx-Jyy*wx)/Jzz 0]
#   return [Jpb;Jw]
# end

# function Api(J::Matrix{T},X::Vector{T})
#   Jxx,Jyy,Jzz = diag(J)
#   px,py,pz = X[1:3]
#   wx,wy,wz = X[4:6]
#   Jpi = [(px*wx)/2+(py*wy)/2+(pz*wz)/2 (px*wy)/2-wz/2-(py*wx)/2 wy/2+(px*wz)/2-(pz*wx)/2 px^2/4-py^2/4-pz^2/4+1/4 pz/2+(px*py)/2 (px*pz)/2-py/2;
#           wz/2-(px*wy)/2+(py*wx)/2 (px*wx)/2+(py*wy)/2+(pz*wz)/2 (py*wz)/2-wx/2-(pz*wy)/2 (px*py)/2-pz/2 -px^2/4+py^2/4-pz^2/4+1/4 px/2+(py*pz)/2;
#           (pz*wx)/2-(px*wz)/2-wy/2 wx/2-(py*wz)/2+(pz*wy)/2 (px*wx)/2+(py*wy)/2+(pz*wz)/2 py/2+(px*pz)/2 (py*pz)/2-px/2 -px^2/4-py^2/4+pz^2/4+1/4]
#   Jw = [0 0 0 0                     (Jyy*wz-Jzz*wz)/Jxx   (Jyy*wy-Jzz*wy)/Jxx;
#         0 0 0 -(Jxx*wz-Jzz*wz)/Jyy  0                     -(Jxx*wx-Jzz*wx)/Jyy;
#         0 0 0 (Jxx*wy-Jyy*wy)/Jzz   (Jxx*wx-Jyy*wx)/Jzz   0]
#   return [Jpi;Jw]
# end


# function B(J::Matrix{T},X::Vector{T})
#   Jxx,Jyy,Jzz = diag(J)
#   Bw = length(X) == 6 ? zeros(T,6,3) : zeros(T,7,3)
#   Bw[end-2:end,:] = diagm([1/Jxx,1/Jyy,1/Jzz])
#   return Bw
# end


# function mrp_jacobian(p,w)
#   # w in inertial frame
#   Ji = zeros(eltype(p), (3,6))
#   p1,p2,p3 = p
#   wi1,wi2,wi3 = w
# 
#   # Jb[1,:] = [(p1*wb1)/2+(p2*wb2)/2+(p3*wb3)/2, wb3/2+(p1*wb2)/2-(p2*wb1)/2, (p1*wb3)/2-wb2/2-(p3*wb1)/2, p1^2/4-p2^2/4-p3^2/4+1/4, (p1*p2)/2-p3/2, p2/2 + (p1*p3)/2]
#   # Jb[2,:] = [(p2*wb1)/2-(p1*wb2)/2-wb3/2, (p1*wb1)/2+(p2*wb2)/2+(p3*wb3)/2, wb1/2+(p2*wb3)/2-(p3*wb2)/2, p3/2+(p1*p2)/2, -p1^2/4+p2^2/4-p3^2/4+1/4, (p2*p3)/2 - p1/2]
#   # Jb[3,:] = [wb2/2-(p1*wb3)/2+(p3*wb1)/2, (p3*wb2)/2-(p2*wb3)/2-wb1/2, (p1*wb1)/2+(p2*wb2)/2+(p3*wb3)/2, (p1*p3)/2-p2/2, p1/2+(p2*p3)/2, -p1^2/4-p2^2/4+p3^2/4+1/4]
#    
#   Ji[1,:] = [(p1*wi1)/2+(p2*wi2)/2+(p3*wi3)/2, (p1*wi2)/2-wi3/2-(p2*wi1)/2, wi2/2+(p1*wi3)/2-(p3*wi1)/2, p1^2/4-p2^2/4-p3^2/4+1/4, p3/2+(p1*p2)/2, (p1*p3)/2-p2/2]
#   Ji[2,:] = [wi3/2-(p1*wi2)/2+(p2*wi1)/2, (p1*wi1)/2+(p2*wi2)/2+(p3*wi3)/2, (p2*wi3)/2-wi1/2-(p3*wi2)/2, (p1*p2)/2-p3/2, -p1^2/4+p2^2/4-p3^2/4+1/4, p1/2+(p2*p3)/2]
#   Ji[3,:] = [(p3*wi1)/2-(p1*wi3)/2-wi2/2, wi1/2-(p2*wi3)/2+(p3*wi2)/2, (p1*wi1)/2+(p2*wi2)/2+(p3*wi3)/2, p2/2+(p1*p3)/2, (p2*p3)/2-p1/2, -p1^2/4-p2^2/4+p3^2/4+1/4]
#   return Ji
# end


# function quat2rpy(q::Vector{T})
#   sinr = 2*(q[4]*q[1] + q[2]*q[3]) 
#   cosr = 1 - 2*(q[1]*q[1] + q[2]*q[2])
#   roll = atan2(sinr, cosr)
# 
#   sinp = 2*(q[4]*q[2] - q[3]*q[1])
#   if abs(sinp) >= 1
#     pitch = sign(sinp)*pi/2
#   else
#     pitch = asin(sinp)
#   end
# 
#   siny = 2*(q[4]*q[3] - q[1]*q[2])
#   cosy = 1 - 2*(q[2]*q[2] + q[3]*q[3])
#   yaw = atan2(siny, cosy)
# 
#   return [roll; pitch; yaw]
# end


# function plot_attitude_sphere(q_matrix::Array, plot_live::Bool=false, meshing::Int=20)
#   if size(q_matrix)[1] == 4
#     q_matrix = q_matrix'
#   end
# 
#   n_samples = size(q_matrix)[1]
# 
#   vector_matrix = zeros(n_samples, 3)
#   [vector_matrix[idx,:] = quat2angle(q_matrix[idx,:]) * quat2vec(q_matrix[idx,:], false) for idx in range(1, n_samples)]
# 
# 
#   # sphere of radius pi
#   phi = linspace(0, 2*pi, meshing)
#   theta = linspace(0, pi, meshing)
#   x = pi*cos(phi)*sin(theta)'
#   y = pi*sin(phi)*sin(theta)'
#   z = pi*ones(meshing)*cos(theta)'
# 
#   figure()
#   plot_wireframe(x,y,z, linestyle="-.")
# 
#   for idx in 1:n_samples
#     scatter3D(vector_matrix[idx, 1], vector_matrix[idx, 2], vector_matrix[idx, 3], color="red")
#     plot_live && pause(1/n_samples)
#   end
# end


function plot_vector_rotation(Qs::Matrix{T}, v::Matrix{T}=eye(T,3)) where T
  # Qs: 4xN matrix of quaternions 
  # v: 3xN vectors to be plotted; defaults to body frame

  n_quat = size(Qs,2)
  n_vec = size(v,2)
  slew_matrix = zeros(T,3,n_quat,n_vec)

  for i in 1:n_quat
    q = Qs[:,i]
    for j in 1:n_vec
      slew_matrix[:,i,j] = quat_rotate(q,v[:,j],"rotation")
    end
  end

  # unit sphere
  meshing = 25
  phi = linspace(0,2*pi,meshing)
  theta = linspace(0,pi,meshing)
  x = cos.(phi)*sin.(theta)'
  y = sin.(phi)*sin.(theta)'
  z = ones(meshing)*cos.(theta)'

  cs = ["red","green","blue"]

  fig = PyPlot.figure()
  for i in 1:n_quat
    clf()
    PyPlot.plot_wireframe(x,y,z, linestyle="-.", color="blue")
    for j in 1:n_vec
      v = slew_matrix[:,i,j]
      PyPlot.plot3D([0; v[1]], [0; v[2]], [0, v[3]], color=cs[j%3+1])
    end
    pause(0.0000001)
  end

  return fig, slew_matrix
end

# function plot_q_slew(q_matrix::Array, v::Array=[-1])
#   # q_matrix: Nx4 matrix
#   # v: 3D vector to be rotated by DCM
# 
#   # 4xN to Nx4
#   if size(q_matrix) == (4,4)
#     # check case of 4x4 matrix and transpose if quaternions are column wise
#     if size(q_matrix)[1] == 4 && (norm(q_matrix[1,:]) - 1 > norm(q_matrix[:,1]) - 1)
#       # NOTE: cannot identify if matrix is identity matrix i.e. eye(4)
#       q_matrix = q_matrix'
#     end
#   elseif size(q_matrix)[1] == 4
#     q_matrix = q_matrix'
#   end
# 
#   numel = size(q_matrix)[1]
# 
#   # basis vectors to slew
#   v1 = [1;0;0]
#   v2 = [0;1;0]
#   v3 = [0;0;1]
# 
#   # generate slew matrix for basis vectors
#   slew_matrix = zeros(numel,3,3)
#   for row_idx in 1:numel
#     slew_matrix[row_idx,:, 1] =  quat_rotate(q_matrix[row_idx,:], v1, "rotation")
#     slew_matrix[row_idx,:, 2] =  quat_rotate(q_matrix[row_idx,:], v2, "rotation")
#     slew_matrix[row_idx,:, 3] =  quat_rotate(q_matrix[row_idx,:], v3, "rotation")
#   end
# 
#   vector_matrix = zeros(numel,3)
#   if v != [-1]
#     v = v/norm(v)
#     for row_idx in 1:numel
#       vector_matrix[row_idx,:] =  quat_rotate(q_matrix[row_idx,:], v, "rotation")
#     end
#   end
# 
#   # unit sphere
#   mag = norm(v)
#   meshing = 50
#   phi = linspace(0, 2*pi, meshing)
#   theta = linspace(0, pi, meshing)
#   x = mag*cos(phi)*sin(theta)'
#   y = mag*sin(phi)*sin(theta)'
#   z = mag*ones(meshing)*cos(theta)'
# 
#   fig1 = figure()
#   for idx in 1:size(slew_matrix)[1]
#     clf()
#     plot_wireframe(x,y,z, linestyle="-.", color="white")
# 
#     v_here = slew_matrix[idx,:,1]
#     plot3D([0; v_here[1]], [0; v_here[2]], [0, v_here[3]], color="red")
# 
#     v_here = slew_matrix[idx,:,2]
#     plot3D([0; v_here[1]], [0; v_here[2]], [0, v_here[3]], color="blue")
# 
#     v_here = slew_matrix[idx,:,3]
#     plot3D([0; v_here[1]], [0; v_here[2]], [0, v_here[3]], color="green")
#     pause(0.00001)
# 
# 
#     if v != [-1]
#       plot3D(vector_matrix[1:idx,1], vector_matrix[1:idx,2], vector_matrix[1:idx,3], color="black")
#     end
#   end
# end
