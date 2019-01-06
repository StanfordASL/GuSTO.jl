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
  clearance

  goal_position_EE_delta_error # Required for numerical imprecisions

  # Parameters that can be updated
  f::Vector
  A::Vector
  B::Vector
end
function PandaKin()
  num_joints = 7
  x_dim = 2*num_joints
  u_dim = num_joints 

  clearance = 0.03
  goal_position_EE_delta_error = 0.01 # 1cm of error admissible for EE position

  PandaKin(x_dim,u_dim,num_joints,[],[],[],[],clearance, goal_position_EE_delta_error, [], [], [])
end

function SCPParam(model::PandaKin, fixed_final_time::Bool)
  convergence_threshold = 0.01
  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::PandaKin)
  Delta0 = 10000. 
  omega0 = 1.
  omegamax = 1.0e10
  epsilon = 1.0e-6
  rho0 = 10. 
  rho1 = 100. 
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
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N
  
  X = hcat(linspace(x_init, x_goal, N)...)
  for i in 1:model.num_joints
    for k=1:N
      sum_squares = sqrt.(X[2*i-1,k]*X[2*i-1,k] + X[2*i,k]*X[2*i,k])
      X[2*i-1,k] = X[2*i-1,k] / sum_squares
      X[2*i,k] = X[2*i,k] / sum_squares
    end
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

  # Update state and control box constraints
  model.x_min = -ones(T,model.x_dim) 
  model.x_max = ones(T,model.x_dim) 
  model.u_min = robot.qd_min
  model.u_max = robot.qd_max

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
    B[2*idx-1, 1] = -x[2*idx]
    B[2*idx,   2] =  x[2*idx-1]
  end
end

## Convex state inequality constraints

## Convex control inequality constraints

## Nonconvex state equality constraints
function ncse_unit_norm_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  X[i,k]^2 + X[i+1,k]^2 - 1. 
end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  # TODO(acauligi): forward kinematics get r for the current obstacle shape
  r = rand(3) 

  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  return clearance - dist
end

function ncsi_EE_goal_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = get_EE_position(x_goal, SCPP)
  p_EE      = get_EE_position(X[:,N], SCPP)

  return p_EE[i] - p_EE_goal[i] - model.goal_position_EE_delta_error
end

function ncsi_EE_goal_constraints_negative(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = get_EE_position(x_goal, SCPP)
  p_EE      = get_EE_position(X[:,N], SCPP)

  return -( p_EE[i] - p_EE_goal[i] ) - model.goal_position_EE_delta_error
end

## Nonconvex state equality constraints (convexified)
function ncse_unit_norm_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  grad = 2*Xp[i:i+1,k] 
  Xp[i,k]^2 + Xp[i+1,k]^2 - 1. + grad'*(X[i:i+1,] - Xp[i:i+1,k])
end

## Nonconvex state inequality constraints (convexified)
function ncsi_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

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

function ncsi_EE_goal_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = get_EE_position(x_goal, SCPP)
  p_EE      = get_EE_position(Xp[:,N], SCPP)
  J_p_EE    = get_EE_Jacobian(Xp[:,N], SCPP)

  return p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]) - p_EE_goal[i] - model.goal_position_EE_delta_error
end

function ncsi_EE_goal_constraints_negative_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_PandaKin(traj, traj_prev, SCPP)

  p_EE_goal = get_EE_position(x_goal, SCPP)
  p_EE      = get_EE_position(Xp[:,N], SCPP)
  J_p_EE    = get_EE_Jacobian(Xp[:,N], SCPP)

  return -( p_EE[i] + transpose(J_p_EE[i,:]) * (X[:,N]-Xp[:,N]) - p_EE_goal[i] ) - model.goal_position_EE_delta_error
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
  # TODO(acauligi): need a way to query workspace position for a given collision object along kinematic chain
  return get_EE_position(traj.X[:,k], SCPP)
end

function SCPConstraints(SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
  model = SCPP.PD.model
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
	for k = 1:N, i = 1:x_dim
    push!(SCPC.convex_state_ineq, (csi_max_bound_constraints, k, i))
    push!(SCPC.convex_state_ineq, (csi_min_bound_constraints, k, i))
	end

  ## Convex control equality constraints
  nothing

  ## Convex control inequality constraints
  for k = 1:N-1
		push!(SCPC.convex_control_ineq, (cci_max_bound_constraints, k, 1))
		push!(SCPC.convex_control_ineq, (cci_min_bound_constraints, k, 1))
  end

  ## Nonconvex state equality constraints
  for k = 1:N, i = 1:2:x_dim
    push!(SCPC.nonconvex_state_eq, (ncse_unit_norm_constraints, k, i))
  end

  ## Nonconvex state inequality constraints
  env_ = WS.btenvironment_keepout
  for k = 1:N, i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_ineq, (ncsi_obstacle_avoidance_constraints, k, i))
  end

  for i = 1:3
    push!(SCPC.nonconvex_state_ineq, (ncsi_EE_goal_constraints         , 0, i))
    push!(SCPC.nonconvex_state_ineq, (ncsi_EE_goal_constraints_negative, 0, i))
  end

  ## Nonconvex state equality constraints (convexified)
  for k = 1:N, i = 1:2:x_dim
    push!(SCPC.nonconvex_state_convexified_eq, (ncse_unit_norm_constraints_convexified, k, i))
  end

  ## Nonconvex state inequality constraints (convexified)
  env_ = WS.btenvironment_keepout
  for k = 1:N, i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_obstacle_avoidance_constraints_convexified, k, i))
  end

  for i = 1:3
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_EE_goal_constraints_convexified         , 0, i))
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_EE_goal_constraints_negative_convexified, 0, i))
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

  # Nonlinear Dynamics
  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k]) + Bp[k]*(U[:,k]-Up[:,k])
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

  # Final EE constraint
  p_EE_goal     = get_EE_position(x_goal, SCPP)
  p_EE_curF     = get_EE_position(X[:,N],  SCPP)
  p_EE_prevF    = get_EE_position(Xp[:,N], SCPP)
  J_p_EE_prevF  = get_EE_Jacobian(Xp[:,N], SCPP)

  num += norm( (p_EE_goal-p_EE_curF) - (p_EE_goal-p_EE_prevF - J_p_EE_prevF*(X[:,N]-Xp[:,N])) )
  den += norm(                         (p_EE_goal-p_EE_prevF - J_p_EE_prevF*(X[:,N]-Xp[:,N])) )

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

function verify_collision_free(traj::Trajectory, SCPP::SCPProblem{PandaBot{T}, PandaKin, E}) where {T,E}
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

function get_configuration(X, model::PandaKin)
  q = zeros(model.num_joints)
  for i in 1:model.num_joints
    q[i] = atan2(X[2*i],X[2*i-1])
  end
  return q
end
get_configuration(X,SCPP::SCPProblem) = get_configuration(X,SCPP.PD.model)

function get_EE_position(X, SCPP::SCPProblem)
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  pan = robot.pan

  q = get_configuration(X,SCPP)
  set_configuration!(robot.state,q)
  p_EE = RigidBodyDynamics.Spatial.transform(robot.state, robot.EE_link_point, robot.world_frame)
  p_EE = p_EE.v

  return p_EE
end

function get_EE_Jacobian(X, SCPP::SCPProblem)
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  pd = robot.pan

  J_p_EE = zeros(3,model.x_dim)

  # J_pEE_joint = dr/dq
  q = get_configuration(X,SCPP)
  set_configuration!(robot.state,q)
  J_pEE_joint = point_jacobian(robot.state,robot.EE_path, RigidBodyDynamics.Spatial.transform(robot.state,robot.EE_link_point,robot.world_frame))
  J_pEE_joint = J_pEE_joint.J[:,1:model.num_joints]

  # J_pEE_embedding = dq/dx
  J_pEE_embedding = zeros(model.num_joints,model.x_dim)
  for i in 1:model.num_joints
    sum_squares = sqrt(X[2*i-1]*X[2*i-1] + X[2*i]*X[2*i])   # TODO(acauligi): does this equal 1?
    J_pEE_embedding[i,2*i-1:2*i] = [-X[2*i]; X[2*i-1]]/sum_squares
  end

  # J = dr/dx = dr/dq * dq/dx
  J_pEE = J_pEE_joint * J_pEE_embedding
  return J_p_EE
end
