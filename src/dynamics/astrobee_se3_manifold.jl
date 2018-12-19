export AstrobeeSE3Manifold
export init_traj_straightline, init_traj_geometricplan

mutable struct AstrobeeSE3Manifold <: DynamicsModel
  # state: r v p omega
  x_dim
  u_dim
  clearance

  x_min
  x_max
  u_min
  u_max

  # Parameters that can be updated
  f::Vector
  A::Vector
  B
end

function AstrobeeSE3Manifold()
  x_dim = 13
  u_dim = 6
  clearance = 0.03

  x_max = 1e3*ones(x_dim)
  x_min = -x_max
  u_max = 1e1*ones(u_dim)
  u_min = -u_max

  AstrobeeSE3Manifold(x_dim, u_dim, clearance, x_min, x_max, u_min, u_max, [], [], [])
end

function SCPParam(model::AstrobeeSE3Manifold, fixed_final_time::Bool)
  convergence_threshold = 0.01
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::AstrobeeSE3Manifold)
  Delta0 = 10.
  omega0 = 1.
  omegamax = 1.0e10
  epsilon = 1.0e-6
  rho0 = 0.01
  rho1 = 0.1
  # alpha_succ = 1.2
  # alpha_fail = 0.5
  beta_succ = 2.
  beta_fail = 0.5
  gamma_fail = 5.

  SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail)
end

function SCPParam_Mao(model::AstrobeeSE3Manifold)
  rho = [0.; 0.25; 0.9]
  Delta_u0 = 1000.
  lambda = 1.
  alpha = 2.
  SCPParam_Mao(rho, Delta_u0, lambda, alpha)
end

function SCPParam_TrajOpt(model::AstrobeeSE3Manifold)
  mu0 = 1.
  s0 = 10.
  c = 10.
  tau_plus = 2. 
  tau_minus = 0.5
  k = 5. 
  ftol = 0.01
  xtol = 0.01
  ctol = 0.01
  max_penalty_iteration = 5
  max_convex_iteration = 5
  max_trust_iteration = 5

  SCPParam_TrajOpt(mu0, s0, c, tau_plus, tau_minus, k, ftol, xtol, ctol,max_penalty_iteration,max_convex_iteration,max_trust_iteration)
end

######
# CVX 
######
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold}) where T
  u_dim = SCPP.PD.model.u_dim

  U,N = traj.U, SCPP.N
  Jm = 0
  for k in 1:N-1
    #Jm += norm(U[:,k])^2
    Jm += sum(U[j,k]^2 for j = 1:u_dim)
  end
  Jm += traj.Tf^2
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  cost_true(traj, traj_prev, SCPP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N

  X = hcat(linspace(x_init, stop=x_goal, length=N)...)
  U = zeros(u_dim, N)
  Trajectory(X, U, tf_guess)
end

# TODO(acauligi): Add geometric plan
function int_traj_geometricplan(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  return Trajectory(TOP)  # Placeholder
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  model.f, model.A, model.B = [], [], B_dyn(Xp[:,1],robot,model)
  for k = 1:N-1
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
    push!(model.A, A_dyn(Xp[:,k],robot,model))
  end
end

function update_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f, A = traj_prev.X, traj_prev.U, model.f, model.A

  for k = 1:N-1
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
    update_A!(A[k], Xp[:,k], robot, model)
  end
end

macro constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
  end
end
# # Same as function above, but wihout state-input which are jump variables.
# macro constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   quote
#     # Todo Thomas: put Tf as a function of the current trajectory, i.d. as a SCPV variable
#     X, U, Tf = $(esc(SCPV)).X, $(esc(SCPV)).U, $(esc(traj_prev)).Tf
#     Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
#     robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
#     x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
#     X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
#   end
# end


## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
  if k == N-1
    return Tf*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    return 0.5*(Tf*(fp + get_f(k+1, model)) + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) +
      Tfp*(Ap*(X[:,k+1]-Xp[:,k+1]) + Bp*(U[:,k+1]-Up[:,k+1]))) - (X[:,k+1]-X[:,k])/dh
  end
end
# function dynamics_constraints_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
# 
#   fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
#   if k == N-1
#     return (  Tf*fp[i] + 
#                 Tfp*( sum(Ap[i,j]*(X[j,k]-Xp[j,k]) for j = 1:x_dim) + sum(Bp[i,j]*(U[j,k]-Up[j,k]) for j = 1:u_dim) ) - 
#                 (X[i,k+1]-X[i,k])/dh
#               )
#   else
#     return (  0.5*( Tf*(fp[i] + get_f(k+1, model)[i]) + 
#                     Tfp*( sum(Ap[i,j]*(X[j,k]-Xp[j,k])     for j = 1:x_dim) + sum(Bp[i,j]*(U[j,k]-Up[j,k])     for j = 1:u_dim)) +
#                     Tfp*( sum(Ap[i,j]*(X[j,k+1]-Xp[j,k+1]) for j = 1:x_dim) + sum(Bp[i,j]*(U[j,k+1]-Up[j,k+1]) for j = 1:u_dim))
#                     ) - 
#               (X[i,k+1]-X[i,k])/dh
#               )
#   end
# end




# Get current dynamics structures for a time step
get_f(k::Int, model::AstrobeeSE3Manifold) = model.f[k]
get_A(k::Int, model::AstrobeeSE3Manifold) = model.A[k]
get_B(k::Int, model::AstrobeeSE3Manifold) = model.B

# Generate full dynamics structures for a time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  x_dim = model.x_dim
  f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  r, v, w = x[1:3], x[4:6], x[11:13]
  qx, qy, qz, qw = x[7:10] 
  wx, wy, wz = x[11:13]
  F, M = u[1:3], u[4:6]

  f[1:3] = v
  f[4:6] = 1/robot.mass*F

  # SO(3)
  f[7] = 1/2*(-wx*qy - wy*qz - wz*qw)
  f[8] = 1/2*( wx*qx - wz*qz + wy*qw)
  f[9] = 1/2*( wy*qx + wz*qy - wx*qw)
  f[10] = 1/2*(wz*qx - wy*qy + wx*qz)
  f[11:13] = robot.Jinv*(M - cross(w,robot.J*w))
end

function A_dyn(x::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)
  A[1:6,1:6] = kron([0 1; 0 0], Eye(3))
  update_A!(A, x, robot, model)
  return A
end

function update_A!(A, x::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  Jxx, Jyy, Jzz = diag(robot.J)
  qx, qy, qz, qw = x[7:10] 
  wx, wy, wz = x[11:13]

  # See scripts/symbolic_math.m for derivation of these equations
  A[7,8] = -wx/2
  A[7,9] = -wy/2
  A[7,10] = -wz/2
  A[7,11] = -qy/2
  A[7,12] = -qz/2
  A[7,13] = -qw/2

  A[8,7] = wx/2
  A[8,9] = -wz/2
  A[8,10] = wy/2
  A[8,11] = qx/2
  A[8,12] = qw/2
  A[8,13] = -qz/2

  A[9,7] = wy/2
  A[9,8] = wz/2
  A[9,10] = -wx/2
  A[9,11] = -qw/2
  A[9,12] = qx/2
  A[9,13] = qy/2

  A[10,7] = wz/2
  A[10,8] = -wy/2
  A[10,9] = wx/2
  A[10,11] = qz/2
  A[10,12] = -qy/2
  A[10,13] = qx/2

  A[11,12] = (Jyy-Jzz)*wz/Jxx
  A[11,13] = (Jyy-Jzz)*wy/Jxx
  A[12,11] = -(Jxx-Jzz)*wz/Jyy
  A[12,13] = -(Jxx-Jzz)*wx/Jyy
  A[13,11] = (Jxx-Jyy)*wy/Jzz
  A[13,12] = (Jxx-Jyy)*wx/Jzz
end

function B_dyn(x::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  x_dim, u_dim = model.x_dim, model.u_dim
  B = zeros(x_dim, u_dim)
  B[1:6,1:3] = 1/robot.mass*kron([0;1], Eye(3))
  B[11:13,4:6] = robot.Jinv   # SO(3)
  return B
end

## Convex state inequality constraints
function csi_orientation_sign(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return -X[10,k]
end

function csi_translational_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  #return norm(X[4:6,k]) - robot.hard_limit_vel
  # Changed from the norm function so it works with JuMP
  return sum(X[3+j,k]^2 for j = 1:3)  - robot.hard_limit_vel^2
end

function csi_angular_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  # Changed from the norm function so it works with JuMP
  #return norm(X[11:13,k]) - robot.hard_limit_omega
  return sum(X[10+j,k]^2 for j = 1:3)  - robot.hard_limit_omega^2
end

# # JuMP
# function csi_orientation_sign_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   return -X[10,k]
# end
# function csi_translational_velocity_bound_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   #return norm(X[4:6,k]) - robot.hard_limit_vel
#   return sum(X[3+j,k]^2 for j = 1:3)  - robot.hard_limit_vel^2
# end
# function csi_angular_velocity_bound_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   #return norm(X[11:13,k]) - robot.hard_limit_omega
#   return sum(X[10+j,k]^2 for j = 1:3)  - robot.hard_limit_omega^2
# end

## Convex control inequality constraints
function cci_translational_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  #return 1/robot.mass*norm(U[1:3,k]) - robot.hard_limit_accel
  return 1/robot.mass*sum(U[j,k]^2 for j = 1:3)  - robot.hard_limit_accel^2
end

function cci_angular_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  #return norm(robot.Jinv*U[4:6,k]) - robot.hard_limit_alpha
  return sum(sum(robot.Jinv[j,i]*U[i,k] for i=1:3)^2 for j = 1:3) - robot.hard_limit_alpha^2
end

# # JuMP
# function cci_translational_accel_bound_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   return 1/robot.mass*sum(U[j,k]^2 for j = 1:3)  - robot.hard_limit_accel^2
# end
# function cci_angular_accel_bound_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
#   return sum(sum(robot.Jinv[j,i]*U[i,k] for i=1:3)^2 for j = 1:3) - robot.hard_limit_alpha^2
# end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  r = X[1:3,k]
  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
  return clearance - dist
end

## Nonconvex state inequality constraints (convexified)
function ncsi_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance

  r0 = get_workspace_location(traj_prev, SCPP, k)
  dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r0, env_idx)

  if dist < SCPP.param.obstacle_toggle_distance
    r = get_workspace_location(traj, SCPP, k)

    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    return (clearance - (dist + nhat'*(r-r0)))
  else
    return 0.
  end
end
function ncsi_obstacle_avoidance_constraints_convexified_jump(SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold_jump(SCPV, traj_prev, SCPP)
  
  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance

  r0 = get_workspace_location(traj_prev, SCPP, k)
  dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r0, env_idx)

  if dist < SCPP.param.obstacle_toggle_distance
    #r = get_workspace_location(traj, SCPP, k)

    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    #return (clearance - (dist + nhat'*(r-r0)))
    return (clearance - sum(nhat[i]*X[i,k] for i=1:3) - nhat'*r0)
  else
    return 0.
  end

  return sum(sum(robot.Jinv[j,i]*U[i,k] for i=1:3)^2 for j = 1:3) - robot.hard_limit_alpha^2
end

## State trust region inequality constraints
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)

  #return norm(X[:,k]-Xp[:,k])^2
  return sum((X[j,k]-Xp[j,k])^2 for j = 1:x_dim)
end

## Trust region inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  #return norm(U[:,k]-Up[:,k])
  return sum((U[j,k]-Up[j,k])^2 for j = 1:u_dim)
end

function get_workspace_location(traj, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int=0) where {T,E}
  return traj.X[1:3,k]
end

# TODO(ambyld): Possibly make a custom version of this for each algorithm
function SCPConstraints(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS = SCPP.WS

  SCPC = SCPConstraints()

  ## Dynamics constraints
  for k = 1:N-1
    push!(SCPC.dynamics, (dynamics_constraints, k, 0))
  end

  ## Convex state equality constraints
  # Init and goal (add init first for convenience of getting dual)
  for i = 1:x_dim
    push!(SCPC.convex_state_eq, (cse_init_constraints, 0, i))
  end
  for i = 1:x_dim
    push!(SCPC.convex_state_eq, (cse_goal_constraints, 0, i))
  end

  ## Convex state inequality constraints
  for k = 1:N
    push!(SCPC.convex_state_ineq, (csi_orientation_sign, k, 0))
    push!(SCPC.convex_state_ineq, (csi_translational_velocity_bound, k, 0))
    push!(SCPC.convex_state_ineq, (csi_angular_velocity_bound, k, 0))
  end

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  # for k = 1:N, i = 1:length(env_.convex_env_components)
  #   push!(SCPC.nonconvex_state_ineq, (ncsi_obstacle_avoidance_constraints, k, i))
  # end

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  env_ = WS.btenvironment_keepout
  # for k = 1:N, i = 1:length(env_.convex_env_components)
  #   push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_obstacle_avoidance_constraints_convexified, k, i))
  # end

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  for k = 1:N-1
    push!(SCPC.convex_control_ineq, (cci_translational_accel_bound, k, 0))
    push!(SCPC.convex_control_ineq, (cci_angular_accel_bound, k, 0))
  end

  ## State trust region ineqality constraints
  for k = 1:N
    push!(SCPC.state_trust_region_ineq, (stri_state_trust_region, k, 0))
  end

  ## Constrol trust region inequality constraints
  for k = 1:N-1
    push!(SCPC.control_trust_region_ineq, (ctri_control_trust_region, k, 0))
  end

  return SCPC
end

# TODO(ambyld): Generalize this
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  fp, Ap = model.f, model.A
  num,den = 0, 0 
  env_ = WS.btenvironment_keepout

  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k])
    num += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den += norm(linearized)
  end
  
  clearance = model.clearance 

  # for k in 1:N
  #   r0 = get_workspace_location(traj, SCPP, k)
  #   r = get_workspace_location(traj_prev, SCPP, k)
  #   for (rb_idx,body_point) in enumerate(env_.convex_robot_components)
  #     for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
  #       dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r0,env_idx)
  #       nhat = dist > 0 ?
  #         (xbody-xobs)./norm(xbody-xobs) :
  #         (xobs-xbody)./norm(xobs-xbody) 
  #       linearized = clearance - (dist + nhat'*(r-r0))
  #       
  #       dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)
  #       
  #       num += abs((clearance-dist) - linearized) 
  #       den += abs(linearized) 
  #     end
  #   end
  # end
  return num/den
end

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  fp = model.f
  num, den = 0, 0
  env_ = WS.btenvironment_keepout
  dt = traj.dt

  for k in 1:N-1
    phi_old = norm(fp[k] - (Xp[:,k]-Xp[:,k])/dtp, 1)
    phi_new = norm(f_dyn(X[:,k],U[:,k],robot,model) - (X[:,k]-X[:,k])/dt, 1)
    phi_hat_new = norm(dynamics_constraints(traj,traj_prev,SCPP,k,0), 1)
    num += (phi_old-phi_new)
    den += (phi_old-phi_hat_new) 
  end

  clearance = model.clearance

  for k in 1:N
    r0 = get_workspace_location(traj, SCPP, k)
    r = get_workspace_location(traj_prev, SCPP, k)
    for (rb_idx,body_point) in enumerate(env_.convex_robot_components)
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r0,env_idx)
        phi_old = clearance-dist
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)
        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 
        phi_new = clearance-dist
        phi_hat_new = clearance - (dist + nhat'*(r-r0))

        num += (phi_old-phi_new)
        den += (phi_old-phi_hat_new) 
      end
    end
  end

  return num/den
end

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  num,den = 0, 0

  cost_true_prev = cost_true(traj_prev, traj_prev, SCPP)
  cost_true_new = cost_true(traj, traj, SCPP)
  for k in 1:N-1
    cost_true_prev += norm(dynamics_constraints(traj_prev, traj_prev, SCPP, k, 0),Inf)
    cost_true_new += norm(dynamics_constraints(traj, traj, SCPP, k, 0),Inf)
  end
  cost_linearized_new = cost_true(traj, traj, SCPP)
  return (cost_true_prev-cost_true_new)/(cost_true_prev-cost_linearized_new)
end


##################
# Shooting Method
##################
function get_dual_cvx(prob::Convex.Problem, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, solver) where {T,E}
	if solver == "Mosek"
		return -MathProgBase.getdual(prob.model)[1:SCPP.PD.model.x_dim]
	else
		return []
	end
end

function get_dual_jump(SCPS::SCPSolution, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
	@show -MathProgBase.getconstrduals(SCPS.solver_model.internalModel)[1:SCPP.PD.model.x_dim]
end
function get_status_jump(SCPS::SCPSolution, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  MathProgBase.status(SCPS.solver_model.internalModel)
end

###########
# Dynamics 
###########
function interpolate_traj(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, dt_min=0.1) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  X, U, dt = traj.X, traj.U, traj.dt

  Nstep = ceil(Int, dt/dt_min)
  Nfull = Nstep*(N-1) + 1
  dtfull = dt/Nstep


  Ufull = Matrix(u_dim, Nfull-1)
  Xfull = Matrix(x_dim, Nfull)
  for k in 1:N-1
    istart = Nstep*(k-1)+1
    Ufull[:, istart:Nstep*k] = repmat(U[:,k], 1, Nstep)

    Xfull[:,istart] = X[:,k]
    Ustep = Ufull[:,istart]
    for i = istart:istart+(Nstep-1)
      k1 = f_dyn(Xfull[:,i],Ustep,robot,model)
      x2 = Xfull[:,i] + 0.5*dtfull*k1
      k2 = f_dyn(x2,Ustep,robot,model)
      x3 = Xfull[:,i] + 0.5*dtfull*k2
      k3 = f_dyn(x3,Ustep,robot,model)
      x4 = Xfull[:,i] + dtfull*k3
      k4 = f_dyn(x4,Ustep,robot,model)
      Xfull[:,i+1] = Xfull[:,i] + 1/6*dtfull*(k1 + 2*k2 + 2*k3 + k4)
    end
  end
  Xfull[:,end] = X[:,end]

  return Trajectory(Xfull, Ufull, Tf)
end

function dynamics_constraint_satisfaction(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  WS = SCPP.WS
  X, U = traj.X, traj.U
  Tf, dt = traj.Tf, traj.dt

  J = 0
  for k in 1:N-1
    J += norm((X[:,k+1]-X[:,k])/dt - f_dyn(X[:,k],U[:,k],robot,model),1)
  end
  return J
end

function verify_collision_free(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model, WS = SCPP.PD.model, SCPP.WS
  N = SCPP.N

  rb_idx = 1
  env_ = WS.btenvironment_keepout

  clearance = model.clearance

  for env_idx = 1:length(env_.convex_env_components)
    for k = 1:N
      r = get_workspace_location(traj, SCPP, k)
      dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r, env_idx)
      if dist < 0
        return false, k, dist
      end
    end
  end

  return true, 0, 0.
end



#######
# JuMP
#######
function cost(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  U,N = traj.U, SCPP.N
  dtp = traj_prev.dt
  Jm = 0
  for k in 1:N-1
    Jm += dtp*norm(U[:,k])^2
  end
  Jm += traj.Tf^2
  return Jm
end

function add_variables!(solver_model::Model, SCPV::SCPVariables{JuMP.Variable}, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  robot, model = SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

  # @variable(solver_model, model.x_min[i] <= X[i=1:x_dim,1:N] <= model.x_max[i])
  # @variable(solver_model, model.u_min[i] <= U[i=1:u_dim,1:N] <= model.u_max[i])
  @variable(solver_model, X[1:x_dim,1:N])
  @variable(solver_model, U[1:u_dim,1:N])
  @variable(solver_model, Tf)

  SCPV.X = X
  SCPV.U = U
  SCPV.Tf = Tf
end


function add_objective!(solver_model::Model, SCPV::SCPVariables, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  robot, model = SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

  U = SCPV.U
  N, dt = SCPP.N, SCPP.tf_guess/SCPP.N

  for i = 1:u_dim
    @NLobjective(solver_model, Min, sum(dt*U[i,k]^2 for k = 1:N-1))
  end
end

function add_constraints!(solver_model::Model, SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  X, U = SCPV.X, SCPV.U
  Xp, Up, Tfp, dtp = traj_prev.X, traj_prev.U, traj_prev.Tf, traj_prev.dt
  Tf = copy(Tfp)
  robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
  x_dim, u_dim, N, dh = model.x_dim, model.u_dim, SCPP.N, SCPP.dh
  WS = SCPP.WS
  
  @constraint(solver_model, SCPV.Tf == Tfp )

  # Init (add these first for shooting method!)
  for i = 1:x_dim
    # @constraint(solver_model, X[i,1] == x_init[i])
    #@constraint(solver_model, 0. == cse_init_constraints_jump(SCPV, traj_prev, SCPP, 0, i) )
    @constraint(solver_model, 0. == cse_init_constraints(SCPV, traj_prev, SCPP, 0, i) )
  end
  # Goal
  for i = 1:x_dim
    # @constraint(solver_model, X[i,N] == x_goal[i])
    #@constraint(solver_model, 0. == cse_goal_constraints_jump(SCPV, traj_prev, SCPP, 0, i) )
    @constraint(solver_model, 0. == cse_goal_constraints(SCPV, traj_prev, SCPP, 0, i) )
  end

  # Dynamics
  # for k = 1:N-1
  #   for i = 1:x_dim
  #     @constraint(solver_model, 0. == dynamics_constraints_jump(SCPV, traj_prev, SCPP, k, i) )
  #   end
  # end
  for k = 1:N-1
    @constraint(solver_model, 0. .== dynamics_constraints(SCPV, traj_prev, SCPP, k, 0) )
  end

  ## Convex state inequality constraints
  for k = 1:N
    #@constraint(solver_model, 0. >= csi_orientation_sign_jump(SCPV, traj_prev, SCPP, k, 0) )
    #@constraint(solver_model, 0. >= csi_translational_velocity_bound_jump(SCPV, traj_prev, SCPP, k, 0) )
    #@constraint(solver_model, 0. >= csi_angular_velocity_bound_jump(SCPV, traj_prev, SCPP, k, 0) )
    @constraint(solver_model, 0. >= csi_orientation_sign(SCPV, traj_prev, SCPP, k, 0) )
    @constraint(solver_model, 0. >= csi_translational_velocity_bound(SCPV, traj_prev, SCPP, k, 0) )
    @constraint(solver_model, 0. >= csi_angular_velocity_bound(SCPV, traj_prev, SCPP, k, 0) )
  end

  ## Convex control inequality constraints
  for k = 1:N-1
    #@constraint(solver_model, 0. >= cci_translational_accel_bound_jump(SCPV, traj_prev, SCPP, k, 0) )
    #@constraint(solver_model, 0. >= cci_angular_accel_bound_jump(SCPV, traj_prev, SCPP, k, 0) )
    @constraint(solver_model, 0. >= cci_translational_accel_bound(SCPV, traj_prev, SCPP, k, 0) )
    @constraint(solver_model, 0. >= cci_angular_accel_bound(SCPV, traj_prev, SCPP, k, 0) )
  end

  # ## Nonconvex state inequality constraints (convexified)
  # env_ = WS.btenvironment_keepout
  # for k = 1:N, i = 1:length(env_.convex_env_components)
  #   @constraint(solver_model, 0. >= ncsi_obstacle_avoidance_constraints_convexified_jump(SCPV, traj_prev, SCPP, k, i) )
  # end
  
end

# function add_nonlinear_constraints!(solver_model::Model, SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
#   X, U = SCPV.X, SCPV.U
#   Xp, Up, dtp = traj_prev.X, traj_prev.U, traj_prev.dt
#   robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
#   x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
# 
#   # Dynamics
#   # for k = 1:N-1
#   #   @NLconstraint(solver_model, X[1,k+1] == X[1,k] + dt*model.v*cos(X[3,k]))
#   #   @NLconstraint(solver_model, X[2,k+1] == X[2,k] + dt*model.v*sin(X[3,k]))
#   #   @constraint(solver_model, X[3,k+1] == X[3,k] + dt*model.k*U[k])
#   # end
#   # Dynamics
#   for k = 1:N-1
#     X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
#     fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
#     for i = 1:x_dim
#       if k == N-1
#       @constraint(solver_model, 0 == 
#           Tf*fp[i] + Tfp*(Ap[i,:]*(X[:,k]-Xp[:,k]) + Bp[i,:]*(U[:,k]-Up[:,k])) - (X[i,k+1]-X[i,k])/dh
#           )
#       else
#         @constraint(solver_model, 0 == 
#           0.5*(Tf*(fp + get_f(k+1, model))[i] + Tfp*(Ap[i,:]*(X[:,k]-Xp[:,k]) + Bp[i,:]*(U[:,k]-Up[:,k])) +
#           Tfp*(Ap[i,:]*(X[:,k+1]-Xp[:,k+1]) + Bp[i,:]*(U[:,k+1]-Up[:,k+1]))) - (X[i,k+1]-X[i,k])/dh
#           )
#       end
#     # @constraint(solver_model, X[1,k+1] == X[1,k] + dtp*model.v*(cos(Xp[3,k]) - sin(Xp[3,k])*(X[3,k] - Xp[3,k])))
#     # @constraint(solver_model, X[2,k+1] == X[2,k] + dtp*model.v*(sin(Xp[3,k]) + cos(Xp[3,k])*(X[3,k] - Xp[3,k])))
#     # @constraint(solver_model, X[3,k+1] == X[3,k] + dtp*model.k*U[k])
#     end
#   end
# 
# 
#   # Init and goal
#   for i = 1:x_dim
#     @constraint(solver_model, X[i,1] == x_init[i])
#     @constraint(solver_model, X[i,N] == x_goal[i])
#   end
# end
