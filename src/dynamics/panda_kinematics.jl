
export init_traj_straightline, init_traj_so1

include("forward_kinematics/joint_functions.jl")

mutable struct PandaKin <: DynamicsModel
  # state: [ (cos(q_i), sin(q_i)) x .. x ... ] (for 7 joint angles q_i, i=1...7)
  x_dim
  u_dim
  num_joints
  
	x_min
	x_max
	u_min
	u_max

  self_clearance
  clearance

  p_EE_goal::Vector
  p_EE_goal_delta_error

  # Parameters that can be updated
  f::Vector
  A::Vector
  B::Vector
end
function PandaKin()
  num_joints = 7
  x_dim = 2*num_joints
  u_dim = num_joints 

  self_clearance = 0.01
  clearance = 0.03

  p_EE_goal = zeros(3)
  p_EE_goal_delta_error = 0.02

  PandaKin(x_dim,u_dim,num_joints,[],[],[],[],self_clearance,clearance,
            p_EE_goal,p_EE_goal_delta_error, [], [], [])
end

function SCPParam(model::PandaKin, fixed_final_time::Bool)
  convergence_threshold = 0.05
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::PandaKin)
  Delta0 = 1. 
  omega0 = 10.
  omegamax = 1.0e10
  epsilon = 1.0e-6
  rho0 = 0.5
  rho1 = 0.9
  beta_succ = 2.
  beta_fail = 0.5
  gamma_fail = 5.

  SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail)
end

function SCPParam_TrajOpt(model::PandaKin)
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
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin}) where T
  N = SCPP.N
  X,U = traj.X, traj.U
  Xp,Up = traj_prev.X, traj_prev.U
  Jm = traj.Tf^2
  for k in 1:N-1
    Jm += norm(U[:,k])^2
  end

  for k in 2:N
    Hess = H_pointing_constraint(Xp[:,k])
    eigv,eigvec = eig(Hess)  
    invalid_idx = find(eigv .<= 0)
    eigv[invalid_idx] = 0.
    Hess = eigvec * diagm(eigv) * eigvec'
    quad_form = 0. 

    if isposdef(-Hess)
      # isposdef() returns false for PSD
      # so check if -Hess is PD instead
      warn("Hessian wasn't PSD!")
    elseif typeof(X[:,k]) == Convex.Variable
      quad_form = quadform(X[:,k]-Xp[:,k], Hess)
    else
      quad_form = (X[:,k]-Xp[:,k])' * Hess * (X[:,k]-Xp[:,k])
    end

    Jm += 0.005*(J_pointing_constraint(Xp[:,k])' * (X[:,k]-Xp[:,k]) + quad_form) 
  end
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  cost_true(traj, traj_prev, SCPP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_straightline(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N

  X = hcat(linspace(x_init, x_goal, N)...)
  U = zeros(u_dim, N)
  Trajectory(X, U, tf_guess)
end

function init_traj_so1(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, robot, x_init, x_goal = TOP.PD.model, TOP.PD.robot, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N
  
  X = hcat(linspace(x_init, x_goal, N)...)
  for k=1:N, i in 1:model.num_joints
    X[2*i-1:2*i,k] ./= norm(X[2*i-1:2*i,k])
  end
  
  th_init, th_goal = get_configuration(x_init,model), get_configuration(x_goal,model)
  th_dot = (th_goal-th_init) / (tf_guess/(N-1))
  U = repmat(th_dot,1,N-1)

  Trajectory(X, U, tf_guess)
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  # Update state and control box constraints
  # TODO(acauligi): consider adding quadrant-wise joint angle limits 
  model.x_min = -ones(T,model.x_dim)
  model.x_max =  ones(T,model.x_dim)

  model.u_min = robot.qd_min
  model.u_max = robot.qd_max

  model.f, model.A, model.B = [], [], []
  for k = 1:N-1
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
    push!(model.A, A_dyn(Xp[:,k],Up[:,k],robot,model))
    push!(model.B, B_dyn(Xp[:,k]        ,robot,model))
  end
end

function update_model_params!(SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f, A, B = traj_prev.X, traj_prev.U, model.f, model.A, model.B

  for k = 1:N-1
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
    update_A!(A[k], Xp[:,k], Up[:,k], robot, model)
    update_B!(B[k], Xp[:,k],          robot, model)
  end
end

macro constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  fk, Ak, Bk = get_f(k, model), get_A(k, model), get_B(k, model)
  if k == N-1
    return Tf*fk + Tfp*(Ak*(X[:,k]-Xp[:,k]) + Bk*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    fkp1, Akp1, Bkp1 = get_f(k+1, model), get_A(k+1, model), get_B(k+1, model)
    return 0.5*(Tf*(fk + fkp1) + Tfp*(Ak*(X[:,k]-Xp[:,k]) + Bk*(U[:,k]-Up[:,k])) +
      Tfp*(Akp1*(X[:,k+1]-Xp[:,k+1]) + Bkp1*(U[:,k+1]-Up[:,k+1]))) - (X[:,k+1]-X[:,k])/dh
  end
end

# Get current dynamics structures for a time step
get_f(k::Int, model::PandaKin) = model.f[k]
get_A(k::Int, model::PandaKin) = model.A[k]
get_B(k::Int, model::PandaKin) = model.B[k]

# Generate full dynamics structures for a time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::PandaKin)
  x_dim = model.x_dim
  f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::PandaKin)
  num_joints = model.num_joints 
  for idx in 1:num_joints
    f[2*idx-1]  = -u[idx]*x[2*idx] 
    f[2*idx]    = u[idx]*x[2*idx-1]
  end
end

function A_dyn(x::Vector, u::Vector, robot::Robot, model::PandaKin)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)
  update_A!(A, x, u, robot, model)
  return A
end

function update_A!(A, x::Vector, u::Vector, robot::Robot, model::PandaKin)
  sparse_list = [sparse([1,2], [2,1], [-u[i],u[i]]) for i in 1:length(u)]
  A = full(blkdiag(sparse_list...))
end

function B_dyn(x::Vector, robot::Robot, model::PandaKin)
  x_dim, u_dim = model.x_dim, model.u_dim
  B = zeros(x_dim, u_dim)
  update_B!(B, x, robot, model)
  return B
end

function update_B!(B, x::Vector, robot::Robot, model::PandaKin)
  for idx in 1:model.num_joints
    B[2*idx-1, idx] = -x[2*idx]
    B[2*idx,   idx] =  x[2*idx-1]
  end
end

## Convex state inequality constraints

## Convex control inequality constraints

## Nonconvex state equality constraints
function ncse_unit_norm_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  X[i,k]^2 + X[i+1,k]^2 - 1. 
end

""" 
Nonconvex state inequality constraints for obstacle avoidance
Inputs: 	k - timestep
			j - index of the robot bubble
			i - index of the environment element 
""" 
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E} 
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  rb_idx, env_idx = j, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  p_joint  = get_joint_position(X[:,k],robot,rb_idx)

  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,p_joint,env_idx)

  return clearance - dist
end

""" 
Nonconvex state inequality constraints for self-collision avoidance
Inputs: 	k - timestep
			j - index of the robot bubble
			i - index of the environment element 
""" 
function ncsi_self_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E} 
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  
  idx1, idx2 = j, i  
  env_ = WS.btenvironment_keepout

  self_clearance = model.self_clearance 

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(X[:,k],model))
  
  bubble_1 = env_.convex_robot_components[idx1]
  r1 = get_bubble_position(robot.bubble_array[idx1],robot)
  BulletCollision.set_transformation(bubble_1, r1)
  
  bubble_2 = env_.convex_robot_components[idx2]
  r2 = get_bubble_position(robot.bubble_array[idx2],robot)
  BulletCollision.set_transformation(bubble_2, r2)
  
  dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)

  return self_clearance - dist
end

function ncsi_EE_goal_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal 
  p_EE      = get_EE_position(X[:,N], robot)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal
  p_EE      = get_EE_position(X[:,N], robot)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min - p_EE[i]
end

## Nonconvex state equality constraints (convexified)
function ncse_unit_norm_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  grad = 2*Xp[i:i+1,k] 
  Xp[i,k]^2 + Xp[i+1,k]^2 - 1. + grad'*(X[i:i+1,] - Xp[i:i+1,k])
end

""" 
Nonconvex state inequality constraints for obstacle avoidance (convexified)
Inputs: 	(...)
		k - timestep
			j - index of the robot bubble
			i - index of the environment element 
""" 
function ncsi_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, 
														 k::Int, j::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  rb_idx, env_idx = j, i
  env_ = WS.btenvironment_keepout
  
  Xq,Xq0 = X[:,k], Xp[:,k]

  clearance = model.clearance 

  p_joint   = get_joint_position(Xq0,robot,rb_idx)
  J_p_joint = get_joint_jacobian(Xq0, robot,rb_idx)

  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,p_joint,env_idx)

  # Convexified obstacle avoidance constraint
  if dist < SCPP.param.obstacle_toggle_distance
    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    return (clearance - (dist + transpose(nhat) * J_p_joint * (Xq-Xq0)))
  else
    return 0.
  end
end

""" 
Nonconvex state inequality constraints for self-collision avoidance (convexified)
Inputs: 	(...)
		k - timestep
			j - index of the robot bubble
			i - index of the environment element 
""" 
function ncsi_self_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, j::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  idx1, idx2 = j, i  
  env_ = WS.btenvironment_keepout

  self_clearance = model.self_clearance 

  Xq,Xq0 = X[:,k], Xp[:,k]

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
  
  bubble_1 = env_.convex_robot_components[idx1]
  r1 = get_bubble_position(robot.bubble_array[idx1],robot)
  BulletCollision.set_transformation(bubble_1, r1)
  
  bubble_2 = env_.convex_robot_components[idx2]
  r2 = get_bubble_position(robot.bubble_array[idx2],robot)
  BulletCollision.set_transformation(bubble_2, r2)
  
  dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)

  # Convexified obstacle avoidance constraint
  if dist < SCPP.param.obstacle_toggle_distance
    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    Jr_∂x = get_bubble_jacobian(robot.bubble_array[idx1],robot) * get_jacobian_embedding(Xq0,SCPP)
    return (self_clearance - (dist + nhat'*Jr_∂x*(Xq-Xq0)))
  else
    return 0.
  end
end

function ncsi_EE_goal_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  
  p_EE_goal = model.p_EE_goal 

  p_EE      = get_EE_position(Xp[:,N], robot)
  J_p_EE    = get_EE_jacobian(Xp[:,N], robot)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]) - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal 

  p_EE      = get_EE_position(Xp[:,N], robot)
  J_p_EE    = get_EE_jacobian(Xp[:,N], robot)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min - (p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]))
end

## State trust region inequality constraints
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  return norm(X[:,k]-Xp[:,k])^2
end

## Trust region inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  return norm(U[:,k]-Up[:,k])
end

function get_workspace_location(traj, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int=0) where {T,E}
  nothing 
end

function SCPConstraints(SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  robot, model = SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS = SCPP.WS

  SCPC = SCPConstraints()

  ## Dynamics constraints
  for k = 1:N-1
    push!(SCPC.dynamics, (dynamics_constraints, k, 0, 0))
  end

  ## Convex state equality constraints
  for i = 1:x_dim
    push!(SCPC.convex_state_eq, (cse_init_constraints, 0, 0, i))
  end

  ## Convex state inequality constraints
  nothing

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  nothing

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  for k = 1:N, j = 1:length(env_.convex_robot_components) , i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_ineq, (ncsi_obstacle_avoidance_constraints, k, j, i))
  end

  ## Nonconvex state equality constraints (convexified)
  nothing 

  ## Nonconvex state inequality constraints (convexified)
  for k = 1:N, j = 1:length(env_.convex_robot_components) , i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_obstacle_avoidance_constraints_convexified, k, j, i))
  end

  for i = 1:3
    push!(SCPC.convex_state_goal_ineq, (ncsi_EE_goal_constraints_convexified         , N, 0, i))
    push!(SCPC.convex_state_goal_ineq, (ncsi_EE_goal_constraints_negative_convexified, N, 0, i))
  end

  ## State trust region ineqality constraints
  for k = 1:N
    push!(SCPC.state_trust_region_ineq, (stri_state_trust_region, k, 0))
  end

  return SCPC
end

function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  fp, Ap, Bp = model.f, model.A, model.B
  
  num_dyn,den_dyn = 0, 0 
  num_con,den_con = 0, 0 
  
  env_ = WS.btenvironment_keepout

  # Nonlinear dynamics 
  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k]) + Bp[k]*(U[:,k]-Up[:,k])
    num_dyn += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den_dyn += norm(linearized)
  end

  # EE pointing constraint
  for k in 2:N
    linearized = get_EE_pointing_constraint(Xp[:,k],robot) + J_pointing_constraint(Xp[:,k])' * (X[:,k]-Xp[:,k])
    num_con += (get_EE_pointing_constraint(X[:,k],robot) - linearized)
    den_con += (linearized)
  end

  # Final EE constraint
  p_EE_goal_min = model.p_EE_goal - model.p_EE_goal_delta_error
  p_EE_goal_max = model.p_EE_goal + model.p_EE_goal_delta_error

  p_EE_prev, J_p_EE_prev  = get_EE_position(Xp[:,N],robot), get_EE_jacobian(Xp[:,N],robot) 
  p_EE_cur                = get_EE_position(X[:,N],robot)
  
  linearized = (p_EE_prev + J_p_EE_prev * (X[:,N]-Xp[:,N])) - p_EE_goal_max
  num_con += sum((p_EE_cur-p_EE_goal_max) - linearized)
  den_con += sum(linearized)

  linearized = p_EE_goal_min - ( p_EE_prev + J_p_EE_prev * (X[:,N]-Xp[:,N]))
  num_con += sum((p_EE_goal_min-p_EE_cur) - linearized)
  den_con += sum(linearized)
  
  # Collision avoidance checks
  clearance, self_clearance = model.clearance, model.self_clearance
  for k in 1:N
    Xq,Xq0 = X[:,k], Xp[:,k]
    
    for (rb_idx_1,body_point) in enumerate(env_.convex_robot_components)
      p_bubble_0    = get_joint_position(Xq0,robot,rb_idx_1)
      J_p_0         = get_joint_jacobian(Xq0,robot,rb_idx_1)
      
      p_bubble      = get_joint_position(Xq,robot,rb_idx_1)
      
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist0, xbody0, xobs0 = BulletCollision.distance(env_, rb_idx_1, p_bubble_0, env_idx)
        nhat = dist0 > 0 ?
          (xbody0-xobs0)./norm(xbody0-xobs0) :
          (xobs0-xbody0)./norm(xobs0-xbody0) 
        linearized = clearance - (dist0 + nhat' * J_p_0 * (Xq-Xq0))
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx_1,p_bubble,env_idx)

        num_con += ((clearance-dist) - linearized) 
        den_con += (linearized) 
      end

      # Check for self-collision
      # TODO(acauligi)
    end
  end

  return (num_dyn + abs(num_con)) / (den_dyn + abs(den_con))
end

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  fp, Ap, Bp = model.f, model.A, model.B
  num,den = 0, 0 

  dt = traj.dt
  env_ = WS.btenvironment_keepout

  # Nonlinear dynamics 
  for k in 1:N-1
    phi_old = norm(fp[k] - (Xp[:,k+1]-Xp[:,k])/dtp, 1)
    phi_new = norm(f_dyn(X[:,k],U[:,k],robot,model) - (X[:,k+1]-X[:,k])/dt, 1)
    phi_hat_new = norm(dynamics_constraints(traj,traj_prev,SCPP,k,0,0), 1)
    num += (phi_old-phi_new)
    den += (phi_old-phi_hat_new) 
  end

  # Final EE constraint
  p_EE_min   = model.p_EE_goal - model.p_EE_goal_delta_error
  p_EE_max   = model.p_EE_goal + model.p_EE_goal_delta_error

  p_EE_prev, J_p_EE_prev  = get_EE_position(Xp[:,N],robot), get_EE_jacobian(Xp[:,N],robot) 
  p_EE_cur                = get_EE_position(X[:,N],robot)

  phi_old     = p_EE_prev - p_EE_max 
  phi_new     = p_EE_cur - p_EE_max 
  phi_hat_new = (p_EE_prev + J_p_EE_prev*(X[:,N]-Xp[:,N])) - p_EE_max
  for i in 1:3
    num += (phi_old[i]-phi_new[i])
    den += (phi_old[i]-phi_hat_new[i]) 
  end
  
  phi_old     = p_EE_min - p_EE_prev
  phi_new     = p_EE_min - p_EE_cur
  phi_hat_new = p_EE_min - (p_EE_prev + J_p_EE_prev*(X[:,N]-Xp[:,N]))
  for i in 1:3
    num += (phi_old[i]-phi_new[i])
    den += (phi_old[i]-phi_hat_new[i]) 
  end

  # Collision avoidance checks
  clearance, self_clearance = model.clearance, model.self_clearance
  for k in 1:N
    Xq,Xq0 = X[:,k], Xp[:,k]
    
    for (rb_idx_1,body_point) in enumerate(env_.convex_robot_components)
      p_bubble_0    = get_joint_position(Xq0,robot,rb_idx_1)
      J_p_0         = get_joint_jacobian(Xq0,robot,rb_idx_1)
      
      p_bubble      = get_joint_position(Xq,robot,rb_idx_1)
      
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist0, xbody0, xobs0 = BulletCollision.distance(env_, rb_idx_1, p_bubble_0, env_idx)
        phi_old = clearance - dist0
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx_1,p_bubble,env_idx)
        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 

        phi_hat_new = clearance - (dist + nhat' * J_p_0 * (Xq-Xq0))
        phi_new = clearance-dist

        num += (phi_old-phi_new)
        den += (phi_old-phi_hat_new)
      end
    end
  end

  return num/den
end


##################
# Shooting Method
##################
function get_dual_cvx(prob::Convex.Problem, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, solver) where {T,E}
	if solver == "Mosek"
		return -MathProgBase.getdual(prob.model)[1:SCPP.PD.model.x_dim]
	else
		return []
	end
end

function get_dual_jump(SCPS::SCPSolution, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
	@show -MathProgBase.getconstrduals(SCPS.solver_model.internalModel)[1:SCPP.PD.model.x_dim]
end


###########
# Dynamics 
###########
function interpolate_traj(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, dt_min=0.1) where {T,E}
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
    # Ufull[:, istart:Nstep*k+1] = collect(linspace()

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

function dynamics_constraint_satisfaction(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
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

function verify_joint_limits(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model

  q_min, q_max = robot.q_min, robot.q_max
  for k in 1:N
    q = get_configuration(traj.X[:,k],model)
    any((q.<q_min) .| (q.>q_max)) && return false
  end
  return true
end

function verify_collision_free(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, WS = SCPP.PD.model, SCPP.WS
  X,U = traj.X, traj.U
  N = SCPP.N
  env_ = WS.btenvironment_keepout

  for k in 1:N
    for rb_idx in 1:length(env_.convex_robot_components)
      p_joint  = get_joint_position(X[:,k],robot,rb_idx)
      for env_idx in 1:length(env_.convex_env_components)
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,p_joint,env_idx)
        if dist < 0
          return false
        end
      end
    end
  end
  return true
end

function verify_EE_goal(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  N, robot, model, WS = SCPP.N, SCPP.PD.robot, SCPP.PD.model, SCPP.WS
  X,U = traj.X, traj.U
  env_ = WS.btenvironment_keepout

  p_EE_goal = model.p_EE_goal 
  p_EE_goal_max = p_EE_goal + model.p_EE_goal_delta_error
  p_EE_goal_min = p_EE_goal - model.p_EE_goal_delta_error

  q = get_configuration(traj.X[:,N],model)
  RigidBodyDynamics.set_configuration!(robot.state,q)
  p_EE = RigidBodyDynamics.Spatial.transform(robot.state, robot.EE_link_point, robot.world_frame)
  p_EE = [val for val in p_EE.v]
  
  return all(p_EE - p_EE_goal_max .< 0) && all(p_EE_goal_min - p_EE .< 0)
end

"""
  x1 = cos(θ), x2 = sin(θ) 	( x1²+x2² ) = 1
    => θ = tan⁻¹(x2/x1) 
"""
function get_configuration(X::Vector, model::PandaKin)
  q = zeros(model.num_joints)
  for i in 1:model.num_joints
    x1,x2 = X[2*i-1:2*i]
    q[i] = atan2(x2,x1)
  end
  return q
end

# """ 
# Computation of the jacobian matrix from the base of the robot to the position of the bubble rb_idx with respect to state variables ([x1=cos(joint_angle); x2=sin(joint_angle); ...])
# 
# Chain-rule to compute the jacobian (from the robot base to the bubble position w.r.t. state variables)
#   	∂r/∂x(x0) = ∂r/∂θ(x0) * ∂θ/∂x(x0)
#   
# Jacobians definitions:
#   Jr_∂x	  := 	∂r/∂x(x0)
#   Jr_∂θ 	:= 	∂r/∂θ(x0)
#   Jθ_∂x 	:= 	∂θ/∂x(x0)
# 
# Computation of Jr_∂θ: 
#   refer to [panda_kinematics.jl::get_bubble_jacobian]
# 
# Computation of Jθ_∂x:
#   x1 = cos(θ), x2 = sin(θ) 	( x1²+x2² ) = 1
#     => θ = tan⁻¹(x2/x1) 
#     => ∂θ/∂x1 = -x2/(x1²+x2²) and ∂θ/∂x1 = -x1/(x1²+x2²)
# """
function get_jacobian_embedding(X::Vector,SCPP::SCPProblem)
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Jθ_∂x = zeros(model.num_joints,model.x_dim)
  for i in 1:model.num_joints
    x1,x2 = X[2*i-1:2*i]
    sum_squares = x1^2 + x2^2     # TODO(acauligi): should this := 1?
    Jθ_∂x[i,2*i-1:2*i] = [1/x2; -1/x1]
  end
  return Jθ_∂x
end

function get_joint_jacobian(X::Vector,robot::PandaBot,idx::Int)
  if     idx == 1; return Jp_panda_joint_1(X);
  elseif idx == 2; return Jp_panda_joint_2(X);
  elseif idx == 3; return Jp_panda_joint_3(X);
  elseif idx == 4; return Jp_panda_joint_4(X);
  elseif idx == 5; return Jp_panda_joint_5(X);
  elseif idx == 6; return Jp_panda_joint_6(X);
  elseif idx == 7; return Jp_panda_joint_7(X);
  elseif idx == 8; return Jp_panda_EE(X);
  else; warn("Invalid index passed to get_joint_jacobian")
  end
end

function get_joint_position(X::Vector,robot::PandaBot,idx::Int)
  if     idx == 1; return p_panda_joint_1(X);
  elseif idx == 2; return p_panda_joint_2(X);
  elseif idx == 3; return p_panda_joint_3(X);
  elseif idx == 4; return p_panda_joint_4(X);
  elseif idx == 5; return p_panda_joint_5(X);
  elseif idx == 6; return p_panda_joint_6(X);
  elseif idx == 7; return p_panda_joint_7(X);
  elseif idx == 8; return p_panda_EE(X);
  else; warn("Invalid index passed to get_joint_position")
  end
end

function get_EE_jacobian(X::Vector,robot::PandaBot)
  Jp_panda_EE(X)
end

function get_EE_position(X::Vector,robot::PandaBot)
  p_panda_EE(X)
end

function get_EE_pointing_constraint(X::Vector,robot::PandaBot)
  pointing_constraint(X)
end

function get_EE_J_pointing_constraint(X::Vector,robot::PandaBot)
  J_pointing_constraint(X)
end
