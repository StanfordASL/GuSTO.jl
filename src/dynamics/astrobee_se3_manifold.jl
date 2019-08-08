export AstrobeeSE3Manifold
export init_traj_straightline

mutable struct AstrobeeSE3Manifold <: DynamicsModel
  x_dim   # state: r v p ω
  u_dim   # control F M
  clearance

  x_min
  x_max

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

  AstrobeeSE3Manifold(x_dim, u_dim, clearance, x_min, x_max, [], [], [])
end

function SCPParam(model::AstrobeeSE3Manifold, fixed_final_time::Bool)
  convergence_threshold = 1e-4 #0.001
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::AstrobeeSE3Manifold)
  Δ0 = 1000.
  ω0 = 1.
  ω_max = 1.0e10
  ε = 1.0e-1
  ρ0 = 0.01
  ρ1 = 100.
  β_succ = 2.
  β_fail = 0.5
  γ_fail = 5.

  SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
end

function SCPParam_Mao(model::AstrobeeSE3Manifold)
  ρ = [0.; 0.25; 0.9]
  Δ_u0 = 1000.
  λ = 1.
  α = 2.
  SCPParam_Mao(ρ, Δ_u0, λ, α)
end

function SCPParam_TrajOpt(model::AstrobeeSE3Manifold)
  μ0 = 1.
  s0 = 10.
  c = 10.
  τ_plus = 2.
  τ_minus = 0.5
  k = 5.
  ftol = 0.01
  xtol = 0.01
  ctol = 0.01
  max_penalty_iteration = 5
  max_convex_iteration = 5
  max_trust_iteration = 5

  SCPParam_TrajOpt(μ0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol, max_penalty_iteration, max_convex_iteration, max_trust_iteration)
end

function cost_true(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E} where {T,E}
  u_dim = OAP.PD.model.u_dim
  U, N, dtp = traj.U, OAP.N, traj_prev.dt
  Jm = 0

  # # Forward Euler
  # for k in 2:N
  #   Jm += sum(dt*U[j,k-1]^2 for j = 1:u_dim)
  # end

  # Trapezoidal
  for k in 2:N
    Jm += sum(1/2*dtp*(U[j,k-1]^2 + U[j,k]^2) for j = 1:u_dim)
  end

  # # Simpson's + Trapezoidal + Forward Euler
  # Jm += sum(1/2*(U[j,2]^2 + U[j,1]^2) for j = 1:u_dim)
  # for k in 3:N-1
  #   Jm += sum(1/2*(U[j,k-1]^2 + U[j,k]^2) for j = 1:u_dim)
  # end
  # Jm += sum(U[j,N-1]^2 for j = 1:u_dim)

  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E} where {T,E}
  cost_true(traj, traj_prev, OAP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_nothing(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

  goal_final = get_first_goal_at_time(goal_set, tf_guess)
  x_goal = goal_final.params.point

  X = repmat(0.5(x_init + x_goal),1,N)
  U = zeros(u_dim,N)
  Trajectory(X, U, tf_guess)
end

function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  
  # Set last state to the center of the goal defined latest in time
  t_goal_final = last(goal_set.goals)[1]
  x_goal = zeros(x_dim)
  for goal in values(inclusive(goal_set.goals, searchequalrange(goal_set.goals, t_goal_final)))
    x_goal[goal.ind_coordinates] = center(goal.params)
  end

  X = hcat(range(x_init, stop=x_goal, length=N)...)
  U = zeros(u_dim, N)
  Trajectory(X, U, tf_guess)
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  model.f, model.A, model.B = [], [], B_dyn(Xp[:,1],robot,model)
  for k = 1:N
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
    push!(model.A, A_dyn(Xp[:,k],robot,model))
  end
end

function update_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f, A = traj_prev.X, traj_prev.U, model.f, model.A

  for k = 1:N
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
    update_A!(A[k], Xp[:,k], robot, model)
  end
end

macro scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, goal_set = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.goal_set
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  
  # p previous
  Xkp, Ukp, Xk = X[:,k-1], U[:,k-1], X[:,k]
  fpkp, Apkp, Bpkp, Xpkp, Upkp = get_f(k-1, model), get_A(k-1, model), get_B(k-1, model), Xp[:,k-1], Up[:,k-1]

  # # Forward Euler rule
  # return (Xkp-Xk) + dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp))

  # # Forward Euler rule
  # if k == N
  #   return (Xkp-Xk) + dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp))
  #   # Trapezoidal rule
  # else
  #   Uk = U[:,k]
  #   fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[:,k]
  #   return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp) +
  #                               fpk  + Apk *(Xk-Xpk)   + Bpk *(Uk-Upk))
  # end

  # Just Trapezoidal rule
  Uk = U[:,k]
  fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[:,k]
  return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp) +
                              fpk  + Apk *(Xk-Xpk)   + Bpk *(Uk-Upk))

  # # Forward Euler rule
  # if k == N
  #   return (Xkp-Xk) + dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp))
  # else
  #   Uk = U[:,k]
  #   fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[:,k]
  #   # Trapezoidal rule
  #   if k == 2
  #     return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp) +
  #                                 fpk  + Apk *(Xk-Xpk)   + Bpk *(Uk-Upk))
  #   # Simpson's rule
  #   else
  #     Xkpp, Ukpp = X[:,k-2], U[:,k-2]
  #     fpkpp, Apkpp, Bpkpp, Xpkpp, Upkpp  = get_f(k-2, model), get_A(k-2, model), get_B(k-2, model), Xp[:,k-2], Up[:,k-2]
  #     return (Xkpp-Xk) + 1/6*dtp.*(fpkpp + Apkpp*(Xkpp-Xpkpp) + Bpkpp*(Ukpp-Upkpp) +
  #                               4*(fpkp  + Apkp *(Xkp-Xpkp)   + Bpkp *(Ukp-Upkp)) +
  #                                  fpk   + Apk  *(Xk-Xpk)     + Bpk  *(Uk-Upk))
  #   end
  # end
end

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
  r, v, ω = x[1:3], x[4:6], x[11:13]
  qw, qx, qy, qz = x[7:10]
  ωx, ωy, ωz = x[11:13]
  F, M = u[1:3], u[4:6]

  f[1:3] = v
  f[4:6] = F/robot.mass

  # SO(3)
  f[7]  = 1/2*(-ωx*qx - ωy*qy - ωz*qz)
  f[8]  = 1/2*( ωx*qw - ωz*qy + ωy*qz)
  f[9]  = 1/2*( ωy*qw + ωz*qx - ωx*qz)
  f[10] = 1/2*( ωz*qw - ωy*qx + ωx*qy)
  f[11:13] = robot.Jinv*(M - cross(ω,robot.J*ω))
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
  qw, qx, qy, qz = x[7:10] 
  ωx, ωy, ωz = x[11:13]

  A[7,8] = -ωx/2
  A[7,9] = -ωy/2
  A[7,10] = -ωz/2
  A[7,11] = -qx/2
  A[7,12] = -qy/2
  A[7,13] = -qz/2

  A[8,7] = ωx/2
  A[8,9] = -ωz/2
  A[8,10] = ωy/2
  A[8,11] = qw/2
  A[8,12] = qz/2
  A[8,13] = -qy/2

  A[9,7] = ωy/2
  A[9,8] = ωz/2
  A[9,10] = -ωx/2
  A[9,11] = -qz/2
  A[9,12] = qw/2
  A[9,13] = qx/2

  A[10,7] = ωz/2
  A[10,8] = -ωy/2
  A[10,9] = ωx/2
  A[10,11] = qy/2
  A[10,12] = -qx/2
  A[10,13] = qw/2

  # TODO: Change to account for nondiagonal inertia
  A[11,12] =  (Jyy-Jzz)*ωz/Jxx
  A[11,13] =  (Jyy-Jzz)*ωy/Jxx
  A[12,11] = -(Jxx-Jzz)*ωz/Jyy
  A[12,13] = -(Jxx-Jzz)*ωx/Jyy
  A[13,11] =  (Jxx-Jyy)*ωy/Jzz
  A[13,12] =  (Jxx-Jyy)*ωx/Jzz
end

function B_dyn(x::Vector, robot::Robot, model::AstrobeeSE3Manifold)
  x_dim, u_dim = model.x_dim, model.u_dim
  B = zeros(x_dim, u_dim)
  B[4:6,1:3] = Eye(3)/robot.mass
  B[11:13,4:6] = robot.Jinv   # SO(3)
  return B
end

## Convex state equality constraints
# Manifold constraint, not used in E-SCP (should be linearized)
function cse_quaternion_norm(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  # return sum(X[6+j,k]^2 for j = 1:4) - 1.0
  qp = Xp[7:10,k] 
  return norm(qp) + sum(qp[j]*(X[6+j,k]-Xp[6+j,k])/norm(qp) for j = 1:4) - 1.0 
end

## Convex state inequality constraints
function csi_orientation_sign(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return -X[7,k]
end

function csi_translational_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return sum(X[3+j,k]^2 for j = 1:3)  - robot.hard_limit_vel^2
end

function csi_angular_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return sum(X[10+j,k]^2 for j = 1:3)  - robot.hard_limit_ω^2
end

## Convex control inequality constraints
function cci_translational_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return 1/robot.mass^2*sum(U[j,k]^2 for j = 1:3)  - robot.hard_limit_accel^2
end

function cci_angular_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return sum(sum(robot.Jinv[j,i]*U[3+i,k] for i=1:3)^2 for j = 1:3) - robot.hard_limit_α^2
end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_signed_distance(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  r = X[1:3,k]
  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
  return clearance - dist
end

function ncsi_obstacle_avoidance_potential(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  obstacle_set = SCPP.PD.env.obstacle_set
  obs = obstacle_set[i]
  r = get_workspace_location(traj, SCPP, k)
  r_prev = get_workspace_location(traj_prev, SCPP, k)

  return obstacle_avoid(r, r_prev, obs, SCPP)
end

function ncsi_keepout_zone_avoidance_potential(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  keepout_zones = SCPP.PD.env.keepout_zones
  koz = keepout_zones[i]
  r = get_workspace_location(traj, SCPP, k)
  r_prev = get_workspace_location(traj_prev, SCPP, k)

  x = obstacle_avoid(r, r_prev, koz, SCPP)

  return x
end

function nsci_sphere_obstacle_potential(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  
  obs = SCPP.PD.env.obstacle_set[i]
  r_obs = origin(obs)
  radius_obs = radius(obs)
  keepout_dist = radius_obs + robot.r
  ε = 1e-1  # Smooth gap from sphere edge
  σ = 1.
  # α = (SCPP.param.alg.ε)/(1/(2*π*σ^2)^(3/2)*exp(-radius_obs^2/(2*σ^2)))
  α = 1e3
  r = get_workspace_location(traj, SCPP, k)

  if norm(r - r_obs) >= keepout_dist + ε
    return 0
  else
    dist = norm(r - r_obs)
    nhat = (r - r_obs)/dist
    normal_dist = α*1/(2*π*σ^2)^(3/2)*exp(-dist^2/(2*σ^2))
    a = smooth_max((keepout_dist + ε)^2 - dist^2)
    b = smooth_max(dist^2 - keepout_dist^2)
    η = a/(a+b)
    c = η*normal_dist
    return c
  end
end

function nsci_sphere_obstacle_potential_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  
  obs = SCPP.PD.env.obstacle_set[i]
  r_obs = origin(obs)
  radius_obs = radius(obs)
  keepout_dist = radius_obs + robot.r
  ε = 1e-1  # Smooth gap from sphere edge
  σ = 1.
  # α = (SCPP.param.alg.ε)/(1/(2*π*σ^2)^(3/2)*exp(-radius_obs^2/(2*σ^2)))
  α = 1e3
  r = get_workspace_location(traj, SCPP, k)
  r_prev = get_workspace_location(traj_prev, SCPP, k)

  if norm(r_prev - r_obs) >= radius_obs + ε
    return 0
  else
    dist = norm(r_prev - r_obs)
    normal_dist = α*1/(2*π*σ^2)^(3/2)*exp(-dist^2/(2*σ^2))
    a = smooth_max((keepout_dist + ε)^2 - dist^2)
    b = smooth_max(dist^2 - keepout_dist^2)
    η = a/(a+b)
    c = η*normal_dist

    # Derivatives
    if (keepout_dist + ε)^2 - dist^2 > 0
      ∇a = ((keepout_dist + ε)^2 - dist^2)^(-2)*(-2*(r_prev-r_obs))*exp(-((keepout_dist + ε)^2 - dist^2)^(-1))
    else
      ∇a = zeros(3)
    end

    if dist^2 - keepout_dist^2 > 0
      ∇b = (dist^2 - keepout_dist^2)^(-2)*(2*(r_prev-r_obs))*exp(-(dist^2 - keepout_dist^2)^(-1))
    else
      ∇b = zeros(3)
    end

    ∇η = ∇a/(a + b) - a/(a + b)^2*(∇a + ∇b)
    ∇normal_dist = α*(2*π*σ^2)^(-3/2)*exp(-1/2*dist^2/σ^2)*(-2*(r_prev - r_obs))/(2*σ^2)

    ∇c = η*∇normal_dist + normal_dist*∇η

    return c + ∇c'*(r-r_prev)
  end
end

function smooth_max(x)
  if x > 0
    return exp(-1/x)
  else
    return zero(x)
  end
end

function obstacle_avoid(r, r_prev, obs::HyperRectangle, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  clearance = SCPP.PD.model.clearance
  radius_robot = SCPP.PD.robot.r

  a_min = minimum(obs) - radius_robot - clearance
  a_max = maximum(obs) + radius_robot + clearance
  a_mid = origin(obs)
  lengths = a_max - a_min

  metric = [dot(1. ./ lengths, r_prev-a_mid); -dot(1. ./ lengths, r_prev-a_mid)]

  i = argmin(metric)
  if i[1] <= 3
    return a_max[i] - r[i]
  else
    return r[i-3] - a_min[i-3]
  end
end

## Nonconvex state inequality constraints (convexified)
function ncsi_obstacle_avoidance_signed_distance_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)

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

## State trust region inequality constraints
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return sum((X[j,k]-Xp[j,k])^2 for j = 1:x_dim)
end

## Trust region inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  return sum((U[j,k]-Up[j,k])^2 for j = 1:u_dim)
end

# function csi_max_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
#   X[i,k] - 5.2
# end

# function csi_min_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int, i::Int) where {T,E}
#   X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
#   -0.2 - X[i,k]
# end

function get_workspace_location(traj, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, k::Int) where {T,E}
  return traj.X[1:3,k]
end

# TODO(ambyld): Possibly make a custom version of this for each algorithm
function SCPConstraints(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS, env = SCPP.WS, SCPP.PD.env
  goal_set, tf_guess = SCPP.PD.goal_set, SCPP.tf_guess

  SCPC = SCPConstraints()

  ## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 2:N)

  ## State init constraints
  add_constraint_category!(SCPC.state_init_eq, sie_init_constraints, :scalar, 1, 1:x_dim)

  ## State boundary condition constraints
  for goal in values(inclusive(goal_set.goals, searchsortedfirst(goal_set.goals, tf_guess), searchsortedlast(goal_set.goals, tf_guess)))
    if typeof(goal.params) == PointGoal
      add_constraint_category!(SCPC.convex_state_boundary_condition_eq, csbce_goal_constraints, goal, :array)
    else
      add_constraint_category!(SCPC.convex_state_boundary_condition_ineq, csbci_goal_constraints, goal, :array)
    end
  end

  ## Convex state equality constraints
  # Manifold constraint, not used in E-SCP
  add_constraint_category!(SCPC.convex_state_eq, cse_quaternion_norm, :scalar, 1:N)

  ## Convex state inequality constraints
  add_constraint_category!(SCPC.convex_state_ineq, csi_orientation_sign, :scalar, 1:N)
  add_constraint_category!(SCPC.convex_state_ineq, csi_translational_velocity_bound, :scalar, 1:N)
  add_constraint_category!(SCPC.convex_state_ineq, csi_angular_velocity_bound, :scalar, 1:N)
  # add_constraint_category!(SCPC.convex_state_ineq, csi_max_bound_constraints, :scalar, 1:N, 1:3)
  # add_constraint_category!(SCPC.convex_state_ineq, csi_min_bound_constraints, :scalar, 1:N, 1:3)

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  N_obs = length(WS.btenvironment_keepout.convex_env_components)
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_obstacle_avoidance_signed_distance, :scalar, 1:N, 1:N_obs)
  # N_obs = length(env.obstacle_set)
  # add_constraint_category!(SCPC.nonconvex_state_ineq, nsci_sphere_obstacle_potential, :scalar, 1:N, 1:N_obs)

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_obstacle_avoidance_signed_distance_convexified, :scalar, 1:N, 1:N_obs)
  
  # N_obs = length(env.obstacle_set)
  # if N_obs > 0
  #   add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_obstacle_avoidance_potential, :scalar, 1:N, 1:N_obs)
  # end
  # N_koz = length(env.keepout_zones)
  # add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_keepout_zone_avoidance_potential, :scalar, 1:N, 1:N_koz)

  # N_obs = length(env.obstacle_set)
  # add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, nsci_sphere_obstacle_potential_convexified, :scalar, 1:N, 1:N_obs)

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints 
  add_constraint_category!(SCPC.convex_control_ineq, cci_translational_accel_bound, :scalar, 1:N-1)
  add_constraint_category!(SCPC.convex_control_ineq, cci_angular_accel_bound, :scalar, 1:N-1)

  ## State trust region ineqality constraints
  # add_constraint_category!(SCPC.state_trust_region_ineq, stri_state_trust_region, :scalar, 1:N)

  ## Constrol trust region inequality constraints
  # add_constraint_category!(SCPC.control_trust_region_ineq, ctri_control_trust_region, :scalar, 1:N-1)

  return SCPC
end

# TODO(ambyld): Generalize this
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  fp, Ap = model.f, model.A
  num,den = 0, 0 
  env_ = WS.btenvironment_keepout

  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k])
    num += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den += norm(linearized)
  end
  
  clearance = model.clearance 

  for k in 1:N
    r0 = get_workspace_location(traj, SCPP, k)
    r = get_workspace_location(traj_prev, SCPP, k)
    for (rb_idx,body_point) in enumerate(env_.convex_robot_components)
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r0,env_idx)
        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 
        linearized = clearance - (dist + nhat'*(r-r0))
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

        num += abs((clearance-dist) - linearized) 
        den += abs(linearized) 
      end
    end
  end
  return num/den

  # num_other, den_other = 0.,0.

  # for k in 1:N
  #   r = get_workspace_location(traj, SCPP, k)
  #   r_prev = get_workspace_location(traj_prev, SCPP, k)
  #   for obs in obstacle_set
  #     r_obs = origin(obs)
  #     radius_obs = radius(obs)
  #     ε = 1e-1  # Smooth gap from sphere edge
  #     σ = 1.
  #     # α = (SCPP.param.alg.ε)/(1/(2*π*σ^2)^(3/2)*exp(-radius_obs^2/(2*σ^2)))
  #     α = 1e3

  #     if norm(r_prev - r_obs) >= radius_obs + ε
  #       c = 0
  #       linearized = 0
  #     else
  #       dist = norm(r_prev - r_obs)
  #       normal_dist = α*1/(2*π*σ^2)^(3/2)*exp(-dist^2/(2*σ^2))
  #       a = smooth_max((radius_obs + ε)^2 - dist^2)
  #       b = smooth_max(dist^2 - radius_obs^2)
  #       η = a/(a+b)
  #       c = η*normal_dist

  #       # Derivatives
  #       if (radius_obs + ε)^2 - dist^2 > 0
  #         ∇a = ((radius_obs + ε)^2 - dist^2)^(-2)*(-2*(r_prev-r_obs))*exp(-((radius_obs + ε)^2 - dist^2)^(-1))
  #       else
  #         ∇a = zeros(3)
  #       end

  #       if dist^2 - radius_obs^2 > 0
  #         ∇b = (dist^2 - radius_obs^2)^(-2)*(2*(r_prev-r_obs))*exp(-(dist^2 - radius_obs^2)^(-1))
  #       else
  #         ∇b = zeros(3)
  #       end

  #       ∇η = ∇a/(a + b) - a/(a + b)^2*(∇a + ∇b)
  #       ∇normal_dist = α*(2*π*σ^2)^(-3/2)*exp(-1/2*dist^2/σ^2)*(-2*(r_prev - r_obs))/(2*σ^2)

  #       ∇c = η*∇normal_dist + normal_dist*∇η

  #       linearized = c + ∇c'*(r-r_prev)
  #     end
  #     num_other += c - linearized
  #     den_other += linearized
  #   end
  # end

  # return (num_dyn + abs(num_other))/(den_dyn + abs(den_other)) # TODO: change this in other models
end

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
  fp = model.f
  num, den = 0, 0
  dt = traj.dt  

  for k in 1:N-1
    ϕ_old = norm(fp[k] - (Xp[:,k]-Xp[:,k])/dtp, 1)
    ϕ_new = norm(f_dyn(X[:,k],U[:,k],robot,model) - (X[:,k]-X[:,k])/dt, 1)
    ϕ_hat_new = norm(dynamics_constraints(traj,traj_prev,SCPP,k,0), 1)
    num += (ϕ_old-ϕ_new)
    den += (ϕ_old-ϕ_hat_new) 
  end

  clearance = model.clearance

  for k in 1:N
    r0 = get_workspace_location(traj, SCPP, k)
    r = get_workspace_location(traj_prev, SCPP, k)
    for (rb_idx,body_point) in enumerate(env_.convex_robot_components)
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r0,env_idx)
        ϕ_old = clearance-dist
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)
        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 
        ϕ_new = clearance-dist
        ϕ_hat_new = clearance - (dist + nhat'*(r-r0))

        num += (ϕ_old-ϕ_new)
        den += (ϕ_old-ϕ_hat_new) 
      end
    end
  end

  return num/den

  # obstacle_set = SCPP.PD.env.obstacle_set

  # for k in 1:N
  #   r = get_workspace_location(traj, SCPP, k)
  #   r_prev = get_workspace_location(traj_prev, SCPP, k)

  #   for obs in obstacle_set
  #     r_obs = origin(obs)
  #     radius_obs = radius(obs)
  #     ε = 1e-1  # Smooth gap from sphere edge
  #     σ = 1.
  #     # α = (SCPP.param.alg.ε)/(1/(2*π*σ^2)^(3/2)*exp(-radius_obs^2/(2*σ^2)))
  #     α = 1e3

  #     if norm(r_prev - r_obs) >= radius_obs + ε
  #       c = 0
  #       linearized = 0
  #     else
  #       dist = norm(r_prev - r_obs)
  #       normal_dist = α*1/(2*π*σ^2)^(3/2)*exp(-dist^2/(2*σ^2))
  #       a = smooth_max((radius_obs + ε)^2 - dist^2)
  #       b = smooth_max(dist^2 - radius_obs^2)
  #       η = a/(a+b)
  #       c = η*normal_dist

  #       # Derivatives
  #       if (radius_obs + ε)^2 - dist^2 > 0
  #         ∇a = ((radius_obs + ε)^2 - dist^2)^(-2)*(-2*(r_prev-r_obs))*exp(-((radius_obs + ε)^2 - dist^2)^(-1))
  #       else
  #         ∇a = zeros(3)
  #       end

  #       if dist^2 - radius_obs^2 > 0
  #         ∇b = (dist^2 - radius_obs^2)^(-2)*(2*(r_prev-r_obs))*exp(-(dist^2 - radius_obs^2)^(-1))
  #       else
  #         ∇b = zeros(3)
  #       end

  #       ∇η = ∇a/(a + b) - a/(a + b)^2*(∇a + ∇b)
  #       ∇normal_dist = α*(2*π*σ^2)^(-3/2)*exp(-1/2*dist^2/σ^2)*(-2*(r_prev - r_obs))/(2*σ^2)

  #       ∇c = η*∇normal_dist + normal_dist*∇η

  #       ϕ_new = c
  #       ϕ_hat_new = c + ∇c'*(r-r_prev)
  #     end
  #     num += (ϕ_old-ϕ_new)
  #     den += (ϕ_old-ϕ_hat_new) 
  #   end
  # end

  # return num/den
end

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3Manifold(traj, traj_prev, SCPP)
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
function get_dual_jump(SCPC::SCPConstraints, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  init_constraint_category = SCPC.state_init_eq[:sie_init_constraints][1]
  -JuMP.dual.([init_constraint_category.con_reference[(i,)] for i = init_constraint_category.ind_other[1]])
end

macro shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  quote
    r, v, ω = $(esc(x))[1:3], $(esc(x))[4:6], $(esc(x))[11:13]
    qw, qx, qy, qz = $(esc(x))[7:10]
    ωx, ωy, ωz = $(esc(x))[11:13]

    pr, pv, pω = $(esc(p))[1:3], $(esc(p))[4:6], $(esc(p))[11:13]
    pqw, pqx, pqy, pqz = $(esc(p))[7:10]

    F, M = $(esc(u))[1:3], $(esc(u))[4:6]

    robot, model, WS, x_init, goal_set = $(esc(SP)).PD.robot, $(esc(SP)).PD.model, $(esc(SP)).WS, $(esc(SP)).PD.x_init, $(esc(SP)).PD.goal_set
  
    r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set
  end
end

function dynamics_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  J, Jinv, mass = robot.J, robot.Jinv, robot.mass

  # State variables
  xdot[1:3] += v                            # rdot
  xdot[4:6] += F/mass                       # vdot
  xdot[7]  += 1/2*(-ωx*qx - ωy*qy - ωz*qz)  # qdot
  xdot[8]  += 1/2*( ωx*qw - ωz*qy + ωy*qz)
  xdot[9]  += 1/2*( ωy*qw + ωz*qx - ωx*qz)
  xdot[10] += 1/2*( ωz*qw - ωy*qx + ωx*qy)
  xdot[11:13] += Jinv*(M - cross(ω,J*ω))    # ωdot

  # Dual variables
  # xdot[14:16] += Zeros(3)
  xdot[17:19] += -pr
  xdot[20] += -1/2*( pqx*ωx + pqy*ωy + pqz*ωz)
  xdot[21] += -1/2*(-pqw*ωx + pqy*ωz - pqz*ωy)
  xdot[22] += -1/2*(-pqw*ωy - pqx*ωz + pqz*ωx)
  xdot[23] += -1/2*(-pqw*ωz + pqx*ωy - pqy*ωx)

  # TODO(ambyld): Revise this
  # a1 = [-J[3,1]*ωy + J[2,1]*ωz
  #        2*J[3,1]*ωx + J[3,2]*ωy - J[1,1]*ωz + J[3,3]*ωz
  #       -2*J[2,1]*ωx + J[1,1]*ωy - J[2,2]*ωy - J[2,3]*ωz]
  # a2 = [-J[3,1]*ωx - 2*J[3,2]*ωy + J[2,2]*ωz - J[3,3]*ωz
  #        J[3,2]*ωx - J[1,2]*ωz
  #        J[1,1]*ωx - J[2,2]*ωx + 2*J[1,2]*ωy + J[1,3]*ωz]
  # a3 = [ J[2,1]*ωx + J[2,2]*ωy - J[3,3]*ωy + 2*J[2,3]*ωz
  #       -J[1,1]*ωx + J[3,3]*ωx - J[1,2]*ωy - 2*J[1,3]*ωz
  #       -J[2,3]*ωx + J[1,3]*ωy]

  # xdot[24] += -dot(pω, Jinv*a1) - 1/2*(-pqw*qx + pqx*qw - pqy*qz + pqz*qy)
  # xdot[25] += -dot(pω, Jinv*a2) - 1/2*(-pqw*qy + pqx*qz + pqy*qw - pqz*qx)
  # xdot[26] += -dot(pω, Jinv*a3) - 1/2*(-pqw*qz - pqx*qy + pqy*qx + pqz*qw)

  xdot[24] += -1/2*(-pqw*qx + pqx*qw - pqy*qz + pqz*qy)
  xdot[25] += -1/2*(-pqw*qy + pqx*qz + pqy*qw - pqz*qx)
  xdot[26] += -1/2*(-pqw*qz - pqx*qy + pqy*qx + pqz*qw)
end

function shooting_ode!(xdot, x, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}, t) where {T,E}
  robot, model = SP.PD.robot, SP.PD.model
  
  x, p = x[1:13], x[14:26]
  u = get_control(x, p, SP)
  
  # Add contributions
  fill!(xdot, 0.)
  dynamics_shooting!(xdot, x, p, u, SP)
  # csi_orientation_sign_shooting!(xdot, x, p, u, SP)
  # csi_translational_velocity_bound_shooting!(xdot, x, p, u, SP)
  # csi_angular_velocity_bound_shooting!(xdot, x, p, u, SP)
  # ncsi_obstacle_avoidance_shooting!(xdot, x, p, u, SP)
  # xdot[14:16] = -obstacle_avoidance_penalty_grad(r, SP)
end

function get_control(x, p, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  robot = SP.PD.robot
  pv, pω = p[4:6,:], p[11:13,:]

  F = pv/(2*robot.mass)
  M = robot.Jinv'*pω/2
  [F; M]
end

function csi_orientation_sign_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  qw < 0. ? xdot[7] += -1. : nothing
end

function csi_translational_velocity_bound_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  if sum(v.^2) > robot.hard_limit_vel^2
    xdot[17:19] += 2*v
  end
end

function csi_angular_velocity_bound_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  if sum(ω.^2) > robot.hard_limit_ω^2
    xdot[24:26] += 2*ω
  end
end

function csi_angular_velocity_bound_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)
  if sum(ω.^2) > robot.hard_limit_ω^2
    xdot[24:26] += 2*ω
  end
end

function nsci_sphere_obstacle_potential_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  r,v,ω,qw,qx,qy,qz,ωx,ωy,ωz,pr,pv,pω,pqw,pqx,pqy,pqz,F,M,robot,model,WS,x_init,goal_set = @shooting_shortcut_AstrobeeSE3Manifold(x, p, u, SP)

  for obs in obstacle_set
    r_obs = origin(obs)
    radius_obs = radius(obs)
    keepout_dist = radius_obs + robot.r
    ε = 1e-1  # Smooth gap from sphere edge
    σ = 1.
    # α = (SCPP.param.alg.ε)/(1/(2*π*σ^2)^(3/2)*exp(-radius_obs^2/(2*σ^2)))
    α = 1e3

    if norm(r - r_obs) >= keepout_dist + ε
      c = 0
    else
      dist = norm(r - r_obs)
      normal_dist = α*1/(2*π*σ^2)^(3/2)*exp(-dist^2/(2*σ^2))
      a = smooth_max((keepout_dist + ε)^2 - dist^2)
      b = smooth_max(dist^2 - keepout_dist^2)
      η = a/(a+b)
      c = η*normal_dist

      # Derivatives
      if (keepout_dist + ε)^2 - dist^2 > 0
        ∇a = ((keepout_dist + ε)^2 - dist^2)^(-2)*(-2*(r-r_obs))*exp(-((keepout_dist + ε)^2 - dist^2)^(-1))
      else
        ∇a = zeros(3)
      end

      if dist^2 - keepout_dist^2 > 0
        ∇b = (dist^2 - keepout_dist^2)^(-2)*(2*(r-r_obs))*exp(-(dist^2 - keepout_dist^2)^(-1))
      else
        ∇b = zeros(3)
      end

      ∇η = ∇a/(a + b) - a/(a + b)^2*(∇a + ∇b)
      ∇normal_dist = α*(2*π*σ^2)^(-3/2)*exp(-1/2*dist^2/σ^2)*(-2*(r_prev - r_obs))/(2*σ^2)

      ∇c = η*∇normal_dist + normal_dist*∇η

      xdot[14:16] += ∇c
    end
  end
end

function ncsi_obstacle_avoidance_shooting!(xdot, x, p, u, SP::ShootingProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  robot, model, env = SP.PD.robot, SP.PD.model, SP.PD.env
  keepout_zones, obstacle_set = env.keepout_zones, env.obstacle_set
  clearance, radius_robot = model.clearance, robot.r

  rdot = xdot[14:16]

  # For hyperrectangles that are not rotated
  # for koz in keepout_zones
  #   a_min = minimum(koz) - radius_robot - clearance
  #   a_max = maximum(koz) + radius_robot + clearance
  #   d = a_max - a_min
  #   if (a_min[1] < r_prev[1] < a_max[1]) && (a_min[2] < r_prev[2] < a_max[2]) && (a_min[3] < r_prev[3] < a_max[3])
  #     dist = [r_prev - a_min; a_max - r_prev]
  #     i = argmin(dist)
  #     if i <= 3
  #       grad[i] += 1.
  #     else
  #       grad[i-3] += -1.
  #     end
  #   end
  # end

  N_obs = length(SP.WS.btenvironment_keepout.convex_env_components)
  env_ = SP.WS.btenvironment_keepout

  rb_idx = 1

  for env_idx = 1:N_obs
    dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, x[1:3], env_idx)
    if 0 <= dist < clearance
      nhat = (xbody-xobs)./norm(xobs-xbody)
      rdot += nhat
    elseif dist < 0
      nhat = (xobs-xbody)./norm(xbody-xobs)
      rdot += nhat
    end
  end
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
