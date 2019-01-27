export PandaKin
export init_traj_straightline, init_traj_geometricplan, init_traj_so1

using SparseArrays

include("forward_kinematics/joint_functions.jl")

mutable struct PandaKin <: DynamicsModel
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

"""
  x1 = cos(θ), x2 = sin(θ) 	( x1²+x2² ) = 1
    => θ = tan⁻¹(x2/x1) 
"""
function get_configuration(X::Vector, model::PandaKin)
  q = zeros(model.num_joints)
  for i in 1:model.num_joints
    x1,x2 = X[2*i-1:2*i]
    q[i] = atan(x2,x1)
  end
  return q
end

function SCPParam(model::PandaKin, fixed_final_time::Bool)
  convergence_threshold = 0.05
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::PandaKin)
  Δ0 = 1.
  ω0 = 10.
  ω_max = 1.0e10
  ε = 1.0e-6
  ρ0 = 0.5
  ρ1 = 0.9
  β_succ = 2.
  β_fail = 0.5
  γ_fail = 5.

  SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
end

function cost_true(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{PandaBot{T}, PandaKin, E} where {T,E}
  u_dim = OAP.PD.model.u_dim

  U, N = traj.U, OAP.N
  Jm = traj.Tf^2
  for k in 1:N-1
    Jm += sum(U[j,k]^2 for j = 1:u_dim)
  end
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{PandaBot{T}, PandaKin, E} where {T,E}
  cost_true(traj, traj_prev, OAP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_nothing(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

  X = repmat(0.5(x_init + x_goal),1,N)
  U = zeros(u_dim,N-1)
  Trajectory(X, U, tf_guess)
end

function init_traj_straightline(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N

  X = hcat(range(x_init, stop=x_goal, length=N)...)
  U = zeros(u_dim, N-1)
  Trajectory(X, U, tf_guess)
end

function init_traj_so1(TOP::TrajectoryOptimizationProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
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

macro scp_shortcut_PandaKin(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)
  fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
  if k == N-1
    return Tf.*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    return 0.5*(Tf.*(fp + get_f(k+1, model)) + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) +
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
    f[2*idx-1]  = -x[2*idx]*u[idx]
    f[2*idx]    =  x[2*idx-1]*u[idx]
  end
end

function A_dyn(x::Vector, u::Vector, robot::Robot, model::PandaKin)
  x_dim = model.x_dim
  A = zeros(x_dim, x_dim)
  update_A!(A, x, u, robot, model)
  return A
end

function update_A!(A, x::Vector, u::Vector, robot::Robot, model::PandaKin)
  sparse_list = [SparseArrays.sparse([1,2], [2,1], [-u[i],u[i]]) for i in 1:length(u)]
  A = Matrix(SparseArrays.blockdiag(sparse_list...))
end

function B_dyn(x::Vector, robot::Robot, model::PandaKin)
  x_dim, u_dim = model.x_dim, model.u_dim
  B = zeros(x_dim, u_dim)
  return B
end

function update_B!(B, x::Vector, robot::Robot, model::PandaKin)
  for idx in 1:model.num_joints
    B[2*idx-1, idx] = -x[2*idx]
    B[2*idx,   idx] =  x[2*idx-1]
  end
end

## Convex state inequality constraints
nothing

## Convex control inequality constraints
nothing

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_signed_distance(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int, j::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

  rb_idx, env_idx = j, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  p_joint  = get_joint_position(X[:,k],robot,rb_idx)

  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,p_joint,env_idx)

  return clearance - dist
end

function ncsi_EE_goal_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)
  p_EE_goal = model.p_EE_goal 
  p_EE      = get_EE_position(X[:,N], robot)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal
  p_EE      = get_EE_position(X[:,N], robot)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min - p_EE[i]
end

## Nonconvex state inequality constraints (convexified)
function ncsi_obstacle_avoidance_signed_distance_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int, j::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

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

function ncsi_EE_goal_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal 

  p_EE      = get_EE_position(Xp[:,N], robot)
  J_p_EE    = get_EE_jacobian(Xp[:,N], robot)

  p_EE_goal_max = p_EE_goal[i] + model.p_EE_goal_delta_error
  return p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]) - p_EE_goal_max 
end

function ncsi_EE_goal_constraints_negative_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = model.p_EE_goal 

  p_EE      = get_EE_position(Xp[:,N], robot)
  J_p_EE    = get_EE_jacobian(Xp[:,N], robot)

  p_EE_goal_min = p_EE_goal[i] - model.p_EE_goal_delta_error
  return p_EE_goal_min - (p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]))
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

## State trust region inequality constraints
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)

  return sum((X[j,k]-Xp[j,k])^2 for j = 1:x_dim)
end

## Trust region inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)
  
  return sum((U[j,k]-Up[j,k])^2 for j = 1:u_dim)
end

function get_workspace_location(traj, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int) where {T,E}
  warn("get_workspace_location() not defined for PandaKin!")
  nothing 
end

# TODO(ambyld): Possibly make a custom version of this for each algorithm
function SCPConstraints(SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model = SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N
  WS, env = SCPP.WS, SCPP.PD.env

  SCPC = SCPConstraints()

  ## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 1:N-1)

  ## Convex state equality constraints
  # Init and goal (add init first for convenience of getting dual)
  add_constraint_category!(SCPC.convex_state_eq, cse_init_constraints, :scalar, 0, 1:x_dim)

  ## Convex state inequality constraints
  nothing

  ## Nonconvex state equality constraints
  nothing

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  N_obs, N_rb = length(env_.convex_env_components), length(env_.convex_robot_components)
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_obstacle_avoidance_signed_distance, :scalar, 1:N-1, 1:N_obs, 1:N_rb)
  
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_EE_goal_constraints, :scalar, 1:N-1, 1:3)
  add_constraint_category!(SCPC.nonconvex_state_ineq, ncsi_EE_goal_constraints_negative, :scalar, 1:N-1, 1:3)

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_obstacle_avoidance_signed_distance_convexified, :scalar, 1:N-1, 1:N_obs, 1:N_rb)

  add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_EE_goal_constraints_convexified, :scalar, 1:N-1, 1:3)
  add_constraint_category!(SCPC.nonconvex_state_convexified_ineq, ncsi_EE_goal_constraints_negative_convexified, :scalar, 1:N-1, 1:3)

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints 
  nothing

  ## State trust region ineqality constraints
  add_constraint_category!(SCPC.state_trust_region_ineq, stri_state_trust_region, :scalar, 1:N-1)

  ## Constrol trust region inequality constraints
  add_constraint_category!(SCPC.control_trust_region_ineq, ctri_control_trust_region, :scalar, 1:N-1)

  return SCPC
end

# TODO(ambyld): Generalize this
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @scp_shortcut_PandaKin(traj, traj_prev, SCPP)
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
      J_p_bubble    = get_joint_jacobian(Xq,robot,rb_idx_1)
      
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

      # # Check for self-collision 
      # for rb_idx_2 in 1:rb_idx_1-1
      #   RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
      #   r0_bubble_2 = get_bubble_position(robot.bubble_array[rb_idx_2],robot)

      #   BulletCollision.set_transformation(bubble_1, r0_bubble_1)
      #   bubble_2 = env_.convex_robot_components[rb_idx_2]
      #   BulletCollision.set_transformation(bubble_2, r0_bubble_2)

      #   dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)
      # 
      #   RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq0,model))
      #   Jr_∂x = get_bubble_jacobian(robot.bubble_array[rb_idx_1],robot) * get_jacobian_embedding(Xq0,SCPP)

      #   nhat = dist > 0 ?
      #     (xbody-xobs)./norm(xbody-xobs) :
      #     (xobs-xbody)./norm(xobs-xbody) 
      #   linearized = clearance - (dist + nhat'*Jr_∂x*(Xq-Xq0))
      #   
      #   RigidBodyDynamics.set_configuration!(robot.state, get_configuration(Xq,model))
      #   r_bubble_1 = get_bubble_position(robot.bubble_array[rb_idx_1],robot)
      #   BulletCollision.set_transformation(bubble_1, r_bubble_1)

      #   r_bubble_2 = get_bubble_position(robot.bubble_array[rb_idx_2],robot)
      #   BulletCollision.set_transformation(bubble_2, r_bubble_2)

      #   dist,xbody,xobs = BulletCollision.distance(bubble_1, bubble_2)

      #   num_con += ((clearance-dist) - linearized) 
      #   den_con += (linearized) 
      # end
    end
  end

  return (num_dyn + abs(num_con)) / (den_dyn + abs(den_con))
end

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  warn("trust_region_ratio_trajopt() not defined for PandaKin!")
  return T(1) 
end

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  warn("trust_region_ratio_mao() not defined for PandaKin!")
  return T(1) 
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

function get_dual_jump(SCPC::SCPConstraints, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  JuMP.dual.([SCPC.convex_state_eq[:cse_init_constraints].con_reference[0,(i,)] for i = SCPC.convex_state_eq[:cse_init_constraints].ind_other[1]])
end

macro shooting_shortcut_PandaKin(x, p, u, SP)
  quote
    robot, model, WS, x_init, x_goal = $(esc(SP)).PD.robot, $(esc(SP)).PD.model, $(esc(SP)).WS, $(esc(SP)).PD.x_init, $(esc(SP)).PD.x_goal

    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14= $(esc(x))[1:end]
    px1,px2,px3,px4,px5,px6,px7,px8,px9,px10,px11,px12,px13,px14= $(esc(p))[1:end]
    th_dot = $(esc(u))[1:end]
  
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,
      px1,px2,px3,px4,px5,px6,px7,px8,px9,px10,px11,px12,px13,px14
      th_dot,robot,model,WS,x_init,x_goal
  end
end

function shooting_ode!(Xdot, X, SP::ShootingProblem{PandaBot{T}, PandaKin, E}, t) where {T,E}
  robot, model = SP.PD.robot, SP.PD.model
  
  x, p, u = X[1:model.x_dim], X[model.x_dim+1:end], zeros(model.u_dim)
  get_control!(u, x, p, SP)
  
  # Add contributions
  fill!(xdot, 0.)
  dynamics_shooting!(xdot, x, p, u, SP)
  obstacle_avoidance_shooting!(xdot, x, p, u, SP)
end

""" 
H = - norm(U,2)^2 + p' * xdot = 
  -q1d^2 - q2d^2 + ... - q7d^2 + (p2*q1d*x1-p1*q1d*x2 + p4*q2d*x3-p3*q2d*x4 +
    ... + p14*q7d*x13-p13*q7d*x14)

xd = dH/dp
pd = -dH/dx = [-p2*q1d;p1*q1d; -p4*q2d;p3*q2d; ...; -p14*q7d;p13*q7d]
"""
function dynamics_shooting!(xdot, x, p, u, SP::ShootingProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,
    px1,px2,px3,px4,px5,px6,px7,px8,px9,px10,px11,px12,px13,px14
    th_dot,robot,model,WS,x_init,x_goal = 
    @shooting_shortcut_PandaKin(x, p, u, SP)

  # State variables
  num_joints = model.num_joints 
  for idx in 1:num_joints
    xdot[2*idx-1]  = -x[2*idx]*u[idx]
    xdot[2*idx]    =  x[2*idx-1]*u[idx]
  end

  # Dual variables
  for idx in 1:num_joints
    xdot[model.x_dim+2*idx-1] = -p[2*idx]*u[idx]
    xdot[model.x_dim+2*idx]   =  p[2*idx-1]*u[idx]
  end
end

""" 
H = - norm(U,2)^2 + p' * xdot = 
  -q1d^2 - q2d^2 + ... - q7d^2 + (p2*q1d*x1-p1*q1d*x2 + p4*q2d*x3-p3*q2d*x4 +
    ... + p14*q7d*x13-p13*q7d*x14)

U = arg min H --> dH/dU = 0
  U = 0.5 * [p2*x1-p1*x2; p4*x3-p3*x4; ...; p14*x13-p13*x14]
"""
function get_control!(u, x, p, SP::ShootingProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  u_dim = length(u)
  for idx in 1:u_dim
    u[idx] = 0.5 * (p[2*idx]*x[2*idx-1] - p[2*idx-1]*x[2*idx])
  end
end

function obstacle_avoidance_shooting!(xdot, x, p, u, SP::ShootingProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  warn("obstacle_avoidance_shooting!() not defined for PandaKin!")
  nothing
end

###########
# Dynamics 
###########
function interpolate_traj(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, dt_min=0.1) where {T,E}
  warn("interpolate_traj() not defined yet for PandaKin!")
  return nothing
end

function dynamics_constraint_satisfaction(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  warn("dynamics_constraint_satisfaction() not defined yet for PandaKin!")
  return nothing
end

function verify_collision_free(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  warn("verify_collision_free() not defined yet for PandaKin!")
  return nothing
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

#######
# JuMP
#######
function add_objective!(solver_model::Model, SCPV::SCPVariables, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  robot, model = SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

  U = SCPV.U
  N, dt = SCPP.N, SCPP.tf_guess/SCPP.N

  @NLobjective(solver_model, Min, sum(dt*U[i,k]^2 for i = 1:u_dim, k = 1:N-1))
end
