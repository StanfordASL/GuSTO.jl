export PandaKin
export init_traj_straightline, init_traj_so1

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
function PandaKin(p_EE_goal::Vector)
  num_joints = 7
  x_dim = 2*num_joints
  u_dim = num_joints 

  self_clearance = 0.01
  clearance = 0.03

  p_EE_goal_delta_error = 0.01

  PandaKin(x_dim,u_dim,num_joints,[],[],[],[],self_clearance,
            clearance,p_EE_goal,p_EE_goal_delta_error, [], [], [])
end
PandaKin() = PandaKin(zeros(3))

function SCPParam(model::PandaKin, fixed_final_time::Bool)
  convergence_threshold = 0.05
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::PandaKin)
  Delta0 = 5. 
  omega0 = 1.
  omegamax = 1.0e10
  epsilon = 1.0e-6
  rho0 = 0.5
  rho1 = 0.9
  beta_succ = 2.
  beta_fail = 0.5
  gamma_fail = 5.

  SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail)
end

######
# CVX 
######
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin}) where T
  U,N = traj.U, SCPP.N
  Jm = 0
  for k in 1:N-1
    Jm += norm(U[:,k])^2
  end
  Jm += traj.Tf^2
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
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)
  fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
  if k == N-1
    return Tf*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    return 0.5*(Tf*(fp + get_f(k+1, model)) + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) +
      Tfp*(Ap*(X[:,k+1]-Xp[:,k+1]) + Bp*(U[:,k+1]-Up[:,k+1]))) - (X[:,k+1]-X[:,k])/dh
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

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(X[:,k],model))
  r = get_bubble_position(robot.bubble_array[rb_idx],robot)

  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

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

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(X[:,N],model))
  p_EE      = get_EE_position(robot)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(X[:,N],model))
  p_EE      = get_EE_position(robot)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min -p_EE[i]
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
  clearance = model.clearance

  Xq,Xq0 = X[:,k], Xp[:,k]

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
  r0 = get_bubble_position(robot.bubble_array[rb_idx],robot)
  dist, xbody, xobs = BulletCollision.distance(env_, rb_idx, r0, env_idx)

  # Convexified obstacle avoidance constraint
  if dist < SCPP.param.obstacle_toggle_distance
    # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
    nhat = dist > 0 ?
      (xbody-xobs)./norm(xbody-xobs) :
      (xobs-xbody)./norm(xobs-xbody)

    Jr_∂x = get_bubble_jacobian(robot.bubble_array[rb_idx],robot) * get_jacobian_embedding(Xq0,SCPP)
    return (clearance - (dist + nhat'*Jr_∂x*(Xq-Xq0)))
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

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xp[:,N],model))
  p_EE      = get_EE_position(robot)
  J_p_EE    = get_EE_jacobian(robot) * get_jacobian_embedding(Xp[:,N],SCPP)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]) - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, j::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal 

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xp[:,N],model))
  p_EE      = get_EE_position(robot)
  J_p_EE    = get_EE_jacobian(robot) * get_jacobian_embedding(Xp[:,N],SCPP)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min - ( p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]))
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
    push!(SCPC.dynamics, (dynamics_constraints, k, 0))
  end

  ## Convex state equality constraints
  for i = 1:x_dim
    push!(SCPC.convex_state_eq, (cse_init_constraints, 0, i))
  end

  ## Convex state inequality constraints
  for k =1:N, i=1:x_dim
    push!(SCPC.convex_state_ineq, (csi_max_bound_constraints, k, 0, i))
    push!(SCPC.convex_state_ineq, (csi_min_bound_constraints, k, 0, i))
  end

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  nothing

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  for k = 1:N, j = 1:length(robot.bubble_array), i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_ineq, (ncsi_obstacle_avoidance_constraints, k, j, i))
  end

  for k = 1:N, j = 1:length(robot.bubble_array), i = 1:j-1
    push!(SCPC.nonconvex_state_ineq, (ncsi_self_avoidance_constraints,   k, j, i))
  end

  for i = 1:3
    push!(SCPC.nonconvex_state_ineq, (ncsi_EE_goal_constraints         , N, 0, i))
    push!(SCPC.nonconvex_state_ineq, (ncsi_EE_goal_constraints_negative, N, 0, i))
  end

  ## Nonconvex state equality constraints (convexified)
  nothing 

  ## Nonconvex state inequality constraints (convexified)
  env_ = WS.btenvironment_keepout
  for k = 1:N, j = 1:length(robot.bubble_array), i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_obstacle_avoidance_constraints_convexified, k, j, i))
  end

  for k = 1:N, j = 1:length(robot.bubble_array), i = 1:j-1
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_self_avoidance_constraints_convexified, k, j, i))
  end

  for i = 1:3
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_EE_goal_constraints_convexified         , N, 0, i))
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_EE_goal_constraints_negative_convexified, N, 0, i))
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
  num,den = 0, 0 
  env_ = WS.btenvironment_keepout

  # Nonlinear dynamics 
  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k]) + Bp[k]*(U[:,k]-Up[:,k])
    num += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den += norm(linearized)
  end
  
  # Collision avoidance checks
  clearance, self_clearance = model.clearance, model.self_clearance
  for k in 1:N
    Xq,Xq0 = X[:,k], Xp[:,k]
    
    for (rb_idx_1,body_point) in enumerate(env_.convex_robot_components)
      bubble_1 = env_.convex_robot_components[rb_idx_1]

      RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq,model))
      r_bubble_1 = get_bubble_position(robot.bubble_array[rb_idx_1],robot)

      RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
      r0_bubble_1    = get_bubble_position(robot.bubble_array[rb_idx_1],robot)
      Jr_∂x = get_bubble_jacobian(robot.bubble_array[rb_idx_1],robot) * get_jacobian_embedding(Xq0,SCPP)

      # Check for collision with obstacles
      for (env_idx,convex_env_component) in enumerate(env_.convex_env_components)
        dist, xbody, xobs = BulletCollision.distance(env_, rb_idx_1, r0_bubble_1, env_idx)
        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 
        linearized = clearance - (dist + nhat'*Jr_∂x*(Xq-Xq0))
        
        dist,xbody,xobs = BulletCollision.distance(env_,rb_idx_1,r_bubble_1,env_idx)

        num += abs((clearance-dist) - linearized) 
        den += abs(linearized) 
      end

      # Check for self-collision 
      for rb_idx_2 in 1:rb_idx_1-1
        RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
        r0_bubble_2 = get_bubble_position(robot.bubble_array[rb_idx_2],robot)

        BulletCollision.set_transformation(bubble_1, r0_bubble_1)
        bubble_2 = env_.convex_robot_components[rb_idx_2]
        BulletCollision.set_transformation(bubble_2, r0_bubble_2)

        dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)
      
        RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
        Jr_∂x = get_bubble_jacobian(robot.bubble_array[rb_idx_1],robot) * get_jacobian_embedding(Xq0,SCPP)

        nhat = dist > 0 ?
          (xbody-xobs)./norm(xbody-xobs) :
          (xobs-xbody)./norm(xobs-xbody) 
        linearized = clearance - (dist + nhat'*Jr_∂x*(Xq-Xq0))
        
        RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq,model))
        r_bubble_1 = get_bubble_position(robot.bubble_array[rb_idx_1],robot)
        BulletCollision.set_transformation(bubble_1, r_bubble_1)

        r_bubble_2 = get_bubble_position(robot.bubble_array[rb_idx_2],robot)
        BulletCollision.set_transformation(bubble_2, r_bubble_2)

        dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)

        num += abs((clearance-dist) - linearized) 
        den += abs(linearized) 
      end
    end
  end
  
  # Final EE constraint
  p_EE_goal     = model.p_EE_goal

  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(X[:,N],model))
  p_EE_curF     = get_EE_position(robot)
  
  RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xp[:,N],model))
  p_EE_prevF    = get_EE_position(robot)
  J_p_EE_prevF  = get_EE_jacobian(robot) * get_jacobian_embedding(Xp[:,N],SCPP)

  linearized = (p_EE_goal-p_EE_prevF - J_p_EE_prevF*(X[:,N]-Xp[:,N]))
  num += norm((p_EE_goal-p_EE_curF) - linearized)
  den += norm(linearized)                         
  
  return (num/den)
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
  N = SCPP.N

  env_ = WS.btenvironment_keepout
  clearance, self_clearance = model.clearance, model.self_clearance

  for k = 1:N
    RigidBodyDynamics.set_configuration!(robot.state, get_configuration(traj.X[:,k],model))

    for rb_idx = 1:length(env_.convex_robot_components)
      bubble_1 = env_.convex_robot_components[rb_idx]
      r1 = get_bubble_position(robot.bubble_array[rb_idx],robot)
      BulletCollision.set_transformation(bubble_1, r1)

      # Check for collision with obstacles
      for env_idx = 1:length(env_.convex_env_components)
        obs = env_.convex_env_components[env_idx] 
        dist,xbody,xobs = BulletCollision.distance(bubble_1, obs)
        if dist < 0.  # TODO(acauligi): should this be clearance?
          return false, k, dist
        end
      end

      # Check for self-collision 
      for idx2 in 1:rb_idx-1 
        bubble_2 = env_.convex_robot_components[idx2]
        # r2 = get_bubble_position(robot.bubble_array[idx2],robot)    # TODO(acauligi): I believe pose of bubble_2 would have
        # BulletCollision.set_transformation(bubble_2, r2)            #   already been updated by this point
        
        dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)
        if dist < 0.  # TODO(acauligi): should this be self_clearance?
          return false, k, dist
        end
      end
    end
  end

  return true, 0., 0.
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

""" 
Computation of the jacobian matrix from the base of the robot to the position of the bubble rb_idx with respect to state variables ([x1=cos(joint_angle); x2=sin(joint_angle); ...])

Chain-rule to compute the jacobian (from the robot base to the bubble position w.r.t. state variables)
  	∂r/∂x(x0) = ∂r/∂θ(x0) * ∂θ/∂x(x0)
  
Jacobians definitions:
  Jr_∂x	  := 	∂r/∂x(x0)
  Jr_∂θ 	:= 	∂r/∂θ(x0)
  Jθ_∂x 	:= 	∂θ/∂x(x0)

Computation of Jr_∂θ: 
  refer to [panda_kinematics.jl::get_bubble_jacobian]

Computation of Jθ_∂x:
  x1 = cos(θ), x2 = sin(θ) 	( x1²+x2² ) = 1
    => θ = tan⁻¹(x2/x1) 
    => ∂θ/∂x1 = -x2/(x1²+x2²) and ∂θ/∂x1 = -x1/(x1²+x2²)
"""
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
