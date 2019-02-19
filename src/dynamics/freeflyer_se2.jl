export FreeflyerSE2
export init_traj_straightline, init_traj_geometricplan

mutable struct FreeflyerSE2 <: DynamicsModel
  x_dim   # state: r v p ω
  u_dim
  clearance

  # Parameters that can be updated
  f::Vector
  A
  B
end

function FreeflyerSE2()
  x_dim = 6
  u_dim = 3
  clearance = 0.05
  FreeflyerSE2(x_dim, u_dim, clearance, [], [], [])
end

function SCPParam(model::FreeflyerSE2, fixed_final_time::Bool)
  convergence_threshold = 0.01
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::FreeflyerSE2)
  Δ0 = 3.
  ω0 = 1.
  ω_max = 1.0e10
  ε = 1.0e-6
  ρ0 = 0.1
  ρ1 = 0.3
  β_succ = 2.
  β_fail = 0.5
  γ_fail = 5.

  SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
end

function SCPParam_Mao(model::FreeflyerSE2)
  ρ = [0.;0.25;0.9]
  Δ_u0 = 0.1
  λ = 1.
  α = 2.
  SCPParam_Mao(ρ, Δ_u0, λ, α)
end

function SCPParam_TrajOpt(model::FreeflyerSE2)
  μ0 = 1.
  s0 = 1.
  c = 10. 
  τ_plus = 2. 
  τ_minus = 0.5 
  k = 5. 
  ftol = 0.01
  xtol = 0.1
  ctol = 0.01
  max_penalty_iteration = 5 
  max_convex_iteration = 5
  max_trust_iteration = 5

  SCPParam_TrajOpt(μ0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol, max_penalty_iteration,max_convex_iteration,max_trust_iteration)
end

function cost_true(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Freeflyer{T}, FreeflyerSE2, E} where {T,E}
  u_dim = OAP.PD.model.u_dim
  U, N, dtp = traj.U, OAP.N, traj_prev.dt
  Jm = 0

  # Trapezoidal
  for k in 2:N
    Jm += sum(1/2*dtp*(U[j,k-1]^2 + U[j,k]^2) for j = 1:u_dim)
  end
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Freeflyer{T}, FreeflyerSE2, E} where {T,E}
  cost_true(traj, traj_prev, OAP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

  x_goal = zeros(x_dim)
  for goal in values(inclusive(goal_set.goals, searchsortedfirst(goal_set.goals, tf_guess), searchsortedlast(goal_set.goals, tf_guess)))
    x_goal[goal.ind_coordinates] = center(goal.params)
  end

  X = hcat(range(x_init, stop=x_goal, length=N)...)
  U = zeros(u_dim, N)

  Trajectory(X, U, tf_guess)
end

# TODO(acauligi): Add geometric plan
function int_traj_geometricplan(TOP::TrajectoryOptimizationProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  return Trajectory(TOP)  # Placeholder
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  model.f, model.A, model.B = [], A_dyn(Xp[:,1],robot,model), B_dyn(Xp[:,1],robot,model)
  for k = 1:N
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
  end

  SCPP.PD.robot.Jcollision = []
  for k = 1:N
    push!(SCPP.PD.robot.Jcollision, [Eye(2) zeros(2,1); zeros(1,3)])
  end
  update_model_params!(SCPP, traj_prev)
end

function update_model_params!(SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f = traj_prev.X, traj_prev.U, model.f

  for k = 1:N
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
  end

  xb, Jcollision = SCPP.PD.robot.xb, SCPP.PD.robot.Jcollision
  for k = 1:N
    thetap = traj_prev.X[3,k]
    Jcollision[k][1,3] = -xb[1]*sin(thetap)-xb[2]*cos(thetap)
    Jcollision[k][2,3] =  xb[1]*cos(thetap)-xb[2]*sin(thetap)
  end
end

macro scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, goal_set = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.goal_set
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)

  Xkp, Ukp, Xk = X[:,k-1], U[:,k-1], X[:,k]
  fpkp, Apkp, Bpkp, Xpkp, Upkp = get_f(k-1, model), get_A(k-1, model), get_B(k-1, model), Xp[:,k-1], Up[:,k-1]

  # Just Trapezoidal rule
  Uk = U[:,k]
  fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[:,k]
  return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp*(Ukp-Upkp) +
                              fpk  + Apk *(Xk-Xpk)   + Bpk *(Uk-Upk))

  # return A_dyn_discrete(X[:,k],dtp,robot,model)*X[:,k] + B_dyn_discrete(X[:,k],dtp,robot,model)*U[:,k] - X[:,k+1]
end

# Get current dynamics structures for a time step
get_f(k::Int, model::FreeflyerSE2) = model.f[k]
get_A(k::Int, model::FreeflyerSE2) = model.A
get_B(k::Int, model::FreeflyerSE2) = model.B

# Generate full dynamics structures for a time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::FreeflyerSE2)
  x_dim = model.x_dim
  f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::FreeflyerSE2)
  F, M = u[1:2], u[3]
  f[4:5] = 1/robot.mass_ff*F
  f[6] = robot.J_ff_inv*M
end

function A_dyn(x::Vector, robot::Robot, model::FreeflyerSE2)
  kron([0 1; 0 0], Eye(3))
end

function B_dyn(x::Vector, robot::Robot, model::FreeflyerSE2)
  B = zeros(6,3)
  B[1:6,1:3] = kron([0;1], Eye(3))
  B[4:5,:] = 1/robot.mass_ff * B[4:5,:]
  B[6,:] = robot.J_ff_inv * B[6,:]
  return B
end

# Generate full discrete update version of dynamics matrices for a time step
# TODO(ambyld): Rename these? Initialize these once?
function A_dyn_discrete(x, dt, robot::Robot, model::FreeflyerSE2)
  kron([1 dt; 0 1], Eye(3))
end

function B_dyn_discrete(x, dt, robot::Robot, model::FreeflyerSE2)
  Jzz_inv = robot.J_ff_inv
  B = [0.5*dt^2*Eye(3); dt*Eye(3)]
  B[1:2,:] *= 1/robot.mass_ff  # translational force -> acceleration
  B[4:5,:] *= 1/robot.mass_ff
  B[3,:] *= Jzz_inv      # rotational moment -> acceleration
  B[6,:] *= Jzz_inv 
  return B
end

function csi_translational_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  return sum(X[3+j,k]^2 for j = 1:2)  - robot.hard_limit_vel^2
end

function csi_angular_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  return X[6,k]^2 - robot.hard_limit_ω^2
end

## Convex control inequality constraints
function cci_translational_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  return 1/robot.mass_ff^2*sum(U[j,k]^2 for j = 1:2)  - robot.hard_limit_accel^2
end

function cci_angular_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  Jzz_inv = robot.J_ff_inv
  return (Jzz_inv*U[3,k])^2 - robot.hard_limit_α^2
end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  r = [X[1:2,k];0]
  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
  return clearance - dist
end


## Nonconvex state inequality constraints (convexified)
function ncsi_body_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)

  env_ = WS.btenvironment_keepout
  env_idx = i
  clearance = model.clearance 

  # Base
  rb_idx = 1
  r0 = get_workspace_location(traj_prev, SCPP, k)
  dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r0, env_idx)
  r = get_workspace_location(traj, SCPP, k)

  if dist < SCPP.param.obstacle_toggle_distance
    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    return clearance - (dist + nhat'*(r-r0))
  else
    return 0.
  end
end

function ncsi_arm_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)

  env_ = WS.btenvironment_keepout
  env_idx = i
  clearance = model.clearance 

  r0 = get_workspace_location(traj_prev, SCPP, k)

  # Arm
  rb_idx = 2
  xb = robot.xb
  thetap = Xp[3,k]
  R01p = [cos(thetap) -sin(thetap) 0.;
          sin(thetap)  cos(thetap) 0.;
          0. 0. 1.]
  r1 = R01p*xb + r0
  dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r1, env_idx)

  # return 0.   # toggle off

  if dist < SCPP.param.obstacle_toggle_distance
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    return clearance - (dist + nhat'*robot.Jcollision[k]*(X[1:3,k]-r0))
  else
    return 0.
  end
end

## State trust region inequality constraints
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  return sum((X[j,k]-Xp[j,k])^2 for j = 1:x_dim)
end

## Convex state inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  return sum((U[j,k]-Up[j,k])^2 for j = 1:u_dim)
end

function get_workspace_location(traj, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}, k::Int, i::Int=0) where {T,E}
  return [traj.X[1:2,k]; 0]
end

function SCPConstraints(SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS = SCPP.WS

  SCPC = SCPConstraints()

  ## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 2:N)

  ## Convex state equality constraints
  # Init and goal (add init first for convenience of getting dual)
  add_constraint_category!(SCPC.convex_state_eq, cse_init_constraints, :scalar, 0, 1:x_dim)
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

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  env_ = WS.btenvironment_keepout
  N_obs = length(WS.btenvironment_keepout.convex_env_components)
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_obstacle_avoidance_signed_distance, :scalar, 1:N, 1:N_obs)

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  add_constraint_category!(SCPC.convex_state_ineq, cci_translational_accel_bound, :scalar, 1:N-1)
  add_constraint_category!(SCPC.convex_state_ineq, cci_angular_accel_bound, :scalar, 1:N-1)

  ## State trust region ineqality constraints
  add_constraint_category!(SCPC.state_trust_region_ineq, stri_state_trust_region, :scalar, 1:N)

  ## Constrol trust region inequality constraints
  add_constraint_category!(SCPC.control_trust_region_ineq, ctri_control_trust_region, :scalar, 1:N-1)

  return SCPC
end

function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
  fp, Ap = model.f, model.A
  num,den = 0, 0 
  env_ = WS.btenvironment_keepout

  for k in 1:N-1
    linearized = fp[k] + Ap*(X[:,k]-Xp[:,k])
    num += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den += norm(linearized)
  end
  
  clearance = model.clearance 

  # TODO(ambyld): Update for arm
  for k in 1:N
    r0 = get_workspace_location(traj, SCPP, k)
    r = get_workspace_location(traj_prev, SCPP, k)
    for (rb_idx, body_point) in enumerate(env_.convex_robot_components)
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

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
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

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @scp_shortcut_freeflyerSE2(traj, traj_prev, SCPP)
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

function get_dual_jump(SCPC::SCPConstraints, SCPP::SCPProblem{Freeflyer{T}, FreeflyerSE2, E}) where {T,E}
  init_constraint_category = SCPC.convex_state_eq[:cse_init_constraints][1]
  -JuMP.dual.([init_constraint_category.con_reference[0,(i,)] for i = init_constraint_category.ind_other[1]])
end

##################
# Dynamics 
##################
function A_SE2(rb::Freeflyer)
  kron([0 1; 0 0], Eye(3))
end

function B_SE2(rb::Freeflyer)
  B = zeros(6,3)
  B[1:6,1:3] = kron([0;1], Eye(3))
  B[4:5,:] = 1/rb.mass_ff * B[4:5,:]
  B[6,:] = rb.J_ff_inv * B[6,:]
  return B
end

function Ak_SE2(rb::Freeflyer, dt)
  kron([1 dt; 0 1], Eye(3))
end

function Bk_SE2(rb::Freeflyer, dt)
  Jzz_inv = robot.J_ff_inv
  Bk = [0.5*dt^2*Eye(3); dt*Eye(3)]
  Bk[1:2,:] *= 1/rb.mass_ff  # translational force -> acceleration
  Bk[4:5,:] *= 1/rb.mass_ff
  Bk[3,:] *= Jzz_inv      # rotational moment -> acceleration
  Bk[6,:] *= Jzz_inv 
  return Bk
end


# continuous time dynamics
function f_SE2(X::Vector, U::Vector, rb::Freeflyer)
  F,M = U[1:2], U[3]
  Jzz_inv = robot.J_ff_inv
  xdot = zeros(length(X))
  xdot[1:3] = xdot[4:6]
  xdot[4:5] = 1/rb.mass_ff*F
  xdot[6] = Jzz_inv*M
  return xdot
end
