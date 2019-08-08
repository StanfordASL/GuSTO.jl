export AstrobeeSE3
export init_traj_straightline

mutable struct AstrobeeSE3 <: DynamicsModel
  # state: r v p ω
  x_dim
  u_dim
  clearance

  # Parameters that can be updated
  f::Vector
  A::Vector
  B
end

function AstrobeeSE3()
  x_dim = 12
  u_dim = 6
  clearance = 0.03
  AstrobeeSE3(x_dim, u_dim, clearance, [], [], [])
end

function SCPParam(model::AstrobeeSE3, fixed_final_time::Bool)
  convergence_threshold = 0.01
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::AstrobeeSE3)
  Δ0 = 10.
  ω0 = 1.
  ω_max = 1.0e10
  ε = 1.0e-6
  ρ0 = 0.01
  ρ1 = 0.05
  β_succ = 2.
  β_fail = 0.5
  γ_fail = 5.

  SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
end

function SCPParam_Mao(model::AstrobeeSE3)
  ρ = [0.; 0.25; 0.9]
  Δ_u0 = 1000.
  λ = 1.
  α = 2.
  SCPParam_Mao(ρ, Δ_u0, λ, α)
end

function SCPParam_TrajOpt(model::AstrobeeSE3)
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

######
# CVX 
######
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3}) where T
  U,N = traj.U, SCPP.N
  Jm = 0
  # Trapezoidal
  for k in 2:N
    Jm += sum(1/2*dtp*(U[j,k-1]^2 + U[j,k]^2) for j = 1:u_dim)
  end
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Astrobee3D{T}, AstrobeeSE3, E} where {T,E}
  cost_true(traj, traj_prev, SCPP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_nothing(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
  model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

  goal_final = get_first_goal_at_time(goal_set, tf_guess)
  x_goal = goal_final.params.point

  X = repmat(0.5(x_init + x_goal),1,N)
  U = zeros(u_dim,N)
  Trajectory(X, U, tf_guess)
end

function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
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
function initialize_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  model.f, model.A, model.B = [], [], B_dyn(Xp[:,1],robot,model)
  for k = 1:N
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
    push!(model.A, A_dyn(Xp[:,k],robot,model))
  end
end

function update_model_params!(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f, A = traj_prev.X, traj_prev.U, model.f, model.A

  for k = 1:N
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
    update_A!(A[k], Xp[:,k], robot, model)
  end
end

macro scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, goal_set = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.goal_set
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)

  # p previous
  Xkp, Ukp, Xk = X[:,k-1], U[:,k-1], X[:,k]
  fpkp, Apkp, Bpkp, Xpkp, Upkp = get_f(k-1, model), get_A(k-1, model), get_B(k-1, model), Xp[:,k-1], Up[:,k-1]

  # Just Trapezoidal rule
  Uk = U[:,k]
  fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[:,k]
  return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp) +
                              fpk  + Apk *(Xk-Xpk)   + Bpk *(Uk-Upk))
end

# Get current dynamics structures for a time step
get_f(k::Int, model::AstrobeeSE3) = model.f[k]
get_A(k::Int, model::AstrobeeSE3) = model.A[k]
get_B(k::Int, model::AstrobeeSE3) = model.B

# Generate full dynamics structures for a time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::AstrobeeSE3)
  x_dim = model.x_dim
  f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::AstrobeeSE3)
  r, v, p, w = x[1:3], x[4:6], x[7:9], x[10:12]
  F, M = u[1:3], u[4:6]

  f[1:3] = v
  f[4:6] = 1/robot.mass*F

  # SO(3)
  f[7:9] = mrp_derivative(p,w)
  f[10:12] = robot.Jinv*(M - cross(w,robot.J*w))
end

function A_dyn(x::Vector, robot::Robot, model::AstrobeeSE3)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)
  A[1:6,1:6] = kron([0 1; 0 0], Eye(3))
  update_A!(A, x, robot, model)
  return A
end

function update_A!(A, x::Vector, robot::Robot, model::AstrobeeSE3)
  Jxx, Jyy, Jzz = diag(robot.J)
  px, py, pz = x[7:9] 
  wx, wy, wz = x[10:12]

  # See scripts/symbolic_math.m for derivation of these equations
  A[7,7] = (px*wx)/2+(py*wy)/2+(pz*wz)/2 
  A[7,8] = wz/2+(px*wy)/2-(py*wx)/2
  A[7,9] = (px*wz)/2-wy/2-(pz*wx)/2
  A[7,10] = px^2/4-py^2/4-pz^2/4+1/4
  A[7,11] = (px*py)/2-pz/2
  A[7,12] = py/2+(px*pz)/2

  A[8,7] = (py*wx)/2-(px*wy)/2-wz/2
  A[8,8] = (px*wx)/2+(py*wy)/2+(pz*wz)/2
  A[8,9] = wx/2+(py*wz)/2-(pz*wy)/2
  A[8,10] = pz/2+(px*py)/2
  A[8,11] = -px^2/4+py^2/4-pz^2/4+1/4
  A[8,12] = (py*pz)/2-px/2

  A[9,7] = wy/2-(px*wz)/2+(pz*wx)/2
  A[9,8] = (pz*wy)/2-(py*wz)/2-wx/2
  A[9,9] = (px*wx)/2+(py*wy)/2+(pz*wz)/2
  A[9,10] = (px*pz)/2-py/2
  A[9,11] = px/2+(py*pz)/2
  A[9,12] = -px^2/4-py^2/4+pz^2/4+1/4

  A[10,11] =  (Jyy-Jzz)*wz/Jxx
  A[10,12] =  (Jyy-Jzz)*wy/Jxx
  A[11,10] = -(Jxx-Jzz)*wz/Jyy
  A[11,12] = -(Jxx-Jzz)*wx/Jyy
  A[12,10] =  (Jxx-Jyy)*wy/Jzz
  A[12,11] =  (Jxx-Jyy)*wx/Jzz
end

function B_dyn(x::Vector, robot::Robot, model::AstrobeeSE3)
  x_dim, u_dim = model.x_dim, model.u_dim
  B = zeros(x_dim, u_dim)
  B[1:6,1:3] = 1/robot.mass*kron([0;1], Eye(3))
  B[10:12,4:6] = robot.Jinv   # SO(3)
  return B
end

## Convex state inequality constraints
function csi_translational_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return sum(X[3+j,k]^2 for j = 1:3)  - robot.hard_limit_vel^2
end

function csi_angular_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return sum(X[9+j,k]^2 for j = 1:3)  - robot.hard_limit_ω^2
end

## Convex control inequality constraints
function cci_translational_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return 1/robot.mass^2*sum(U[j,k]^2 for j = 1:3)  - robot.hard_limit_accel^2
end

function cci_angular_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return sum(sum(robot.Jinv[j,i]*U[3+i,k] for i=1:3)^2 for j = 1:3) - robot.hard_limit_α^2
end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  r = X[1:3,k]
  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
  return clearance - dist
end

## Nonconvex state inequality constraints (convexified)
function ncsi_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)

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
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return sum((X[j,k]-Xp[j,k])^2 for j = 1:x_dim)
end

## Trust region inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  return sum((U[j,k]-Up[j,k])^2 for j = 1:u_dim)
end

function get_workspace_location(traj, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, k::Int, i::Int=0) where {T,E}
  return traj.X[1:3,k]
end

# TODO(ambyld): Possibly make a custom version of this for each algorithm
function SCPConstraints(SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
  model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS = SCPP.WS
  goal_set, tf_guess = SCPP.PD.goal_set, SCPP.tf_guess

  SCPC = SCPConstraints()

  ## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 2:N)

  ## State init constraints
  add_constraint_category!(SCPC.state_init_eq, sie_init_constraints, :scalar, 1, 1:x_dim)

  ## State boundary condition constraints
  # add_constraint_category!(SCPC.convex_state_boundary_condition_eq, csbce_goal_constraints, get_first_goal_at_time(goal_set, tf_guess), :array)
  for goal in values(inclusive(goal_set.goals, searchsortedfirst(goal_set.goals, tf_guess), searchsortedlast(goal_set.goals, tf_guess)))
    if typeof(goal.params) == PointGoal
      add_constraint_category!(SCPC.convex_state_boundary_condition_eq, csbce_goal_constraints, goal, :array)
    else
      add_constraint_category!(SCPC.convex_state_boundary_condition_ineq, csbci_goal_constraints, goal, :array)
    end
  end

  ## Convex state inequality constraints
  add_constraint_category!(SCPC.convex_state_ineq, csi_translational_velocity_bound, :scalar, 1:N)
  add_constraint_category!(SCPC.convex_state_ineq, csi_angular_velocity_bound, :scalar, 1:N)

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  N_obs = length(WS.btenvironment_keepout.convex_env_components)
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_obstacle_avoidance_signed_distance, :scalar, 1:N, 1:N_obs)

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_obstacle_avoidance_signed_distance_convexified, :scalar, 1:N, 1:N_obs)

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  add_constraint_category!(SCPC.convex_control_ineq, cci_translational_accel_bound, :scalar, 1:N-1)
  add_constraint_category!(SCPC.convex_control_ineq, cci_angular_accel_bound, :scalar, 1:N-1)

  ## State trust region ineqality constraints
  add_constraint_category!(SCPC.state_trust_region_ineq, stri_state_trust_region, :scalar, 1:N)

  ## Constrol trust region inequality constraints
  add_constraint_category!(SCPC.control_trust_region_ineq, ctri_control_trust_region, :scalar, 1:N-1)

  return SCPC
end

# TODO(ambyld): Generalize this
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
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
end

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
  fp = model.f
  num, den = 0, 0
  env_ = WS.btenvironment_keepout
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
end

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_AstrobeeSE3(traj, traj_prev, SCPP)
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
function get_dual_cvx(prob::Convex.Problem, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, solver) where {T,E}
	if solver == "Mosek"
		return -MathProgBase.getdual(prob.model)[1:SCPP.PD.model.x_dim]
	else
		return []
	end
end

function get_dual_jump(SCPS::SCPSolution, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
	@show -MathProgBase.getconstrduals(SCPS.solver_model.internalModel)[1:SCPP.PD.model.x_dim]
end

###########
# Dynamics 
###########
function interpolate_traj(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}, dt_min=0.1) where {T,E}
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

function dynamics_constraint_satisfaction(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
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

function verify_collision_free(traj::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3, E}) where {T,E}
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
