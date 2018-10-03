export AstrobeeSE2
export init_traj_straightline, init_traj_geometricplan

mutable struct AstrobeeSE2 <: DynamicsModel
  # state: r v p omega
  x_dim
  u_dim
  clearance

  # Parameters that can be updated
  f::Vector
  A
  B
end

function AstrobeeSE2()
  x_dim = 6
  u_dim = 3
  clearance = 0.05
  AstrobeeSE2(x_dim, u_dim, clearance, [], [], [])
end

function SCPParam(model::AstrobeeSE2, fixed_final_time::Bool)
  convergence_threshold = 0.01

  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::AstrobeeSE2)
  Delta0 = 2.
  omega0 = 1.
  omegamax = 1.0e10
  epsilon = 1.0e-6
  rho0 = 0.01
  rho1 = 0.3
  beta_succ = 2.
  beta_fail = 0.5
  gamma_fail = 5.

  SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail)
end

function SCPParam_Mao(model::AstrobeeSE2)
  rho = [0.;0.25;0.9]
  Delta_u0 = 0.1
  lambda = 1.
  alpha = 2.
  SCPParam_Mao(rho, Delta_u0, lambda, alpha)
end

function SCPParam_TrajOpt(model::AstrobeeSE2)
  mu0 = 1.
  s0 = 1.
  c = 10. 
  tau_plus = 2. 
  tau_minus = 0.5 
  k = 5. 
  ftol = 0.01
  xtol = 0.1
  ctol = 0.01
  max_penalty_iteration = 5 
  max_convex_iteration = 5
  max_trust_iteration = 5

  SCPParam_TrajOpt(mu0, s0, c, tau_plus, tau_minus, k, ftol, xtol, ctol,max_penalty_iteration,max_convex_iteration,max_trust_iteration)
end

###############
# Gurobi stuff
###############
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2}) where T
  U, N = traj.U, SCPP.N
  dtp = traj_prev.dt
  Jm = 0
  for k in 1:N-1
    Jm += dtp*norm(U[:,k])^2
  end
  return Jm
end

function cost_true_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  cost_true(traj, traj_prev, SCPP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
  x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess
  N = TOP.N

  X = hcat(linspace(x_init, x_goal, N)...)
  U = zeros(u_dim, N)
  Trajectory(X, U, tf_guess)
end

# TODO(acauligi): Add geometric plan
function int_traj_geometricplan(TOP::TrajectoryOptimizationProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  return Trajectory(TOP)  # Placeholder
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim = model.x_dim, model.u_dim
  Xp, Up = traj_prev.X, traj_prev.U

  model.f, model.A, model.B = [], A_dyn(Xp[:,1],robot,model), B_dyn(Xp[:,1],robot,model)
  for k = 1:N-1
    push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
  end

  SCPP.PD.robot.Jcollision = []
  Jcollision = SCPP.PD.robot.Jcollision
  for k = 1:N
    push!(Jcollision, [eye(2) zeros(2,1); zeros(1,3)])
  end
end

function update_model_params!(SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, traj_prev::Trajectory) where {T,E}
  N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
  Xp, Up, f = traj_prev.X, traj_prev.U, model.f

  for k = 1:N-1
    update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
  end

  xb, Jcollision = SCPP.PD.robot.xb, SCPP.PD.robot.Jcollision
  for k = 1:N
    thetap = traj_prev.X[3,k]
    Jcollision[k][1,3] = -xb[1]*sin(thetap)-xb[2]*cos(thetap)
    Jcollision[k][2,3] =  xb[1]*cos(thetap)-xb[2]*sin(thetap)
  end
end

macro constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  quote
    X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
  end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)

  # Exact dynamics propagation TODO(acauligi): Which version to keep?
  return A_dyn_discrete(X[:,k],dtp,robot,model)*X[:,k] + B_dyn_discrete(X[:,k],dtp,robot,model)*U[:,k] - X[:,k+1]

  fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
  if k == N-1
    return Tf*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    return 0.5*(Tf*(fp + get_f(k+1, model)) + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[:,k]-Up[:,k])) +
      Tfp*(Ap*(X[:,k+1]-Xp[:,k+1]) + Bp*(U[:,k+1]-Up[:,k+1]))) - (X[:,k+1]-X[:,k])/dh
  end
end

# Get current dynamics structures for a time step
get_f(k::Int, model::AstrobeeSE2) = model.f[k]
get_A(k::Int, model::AstrobeeSE2) = model.A
get_B(k::Int, model::AstrobeeSE2) = model.B

# Generate full dynamics structures for a time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::AstrobeeSE2)
  x_dim = model.x_dim
  f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::AstrobeeSE2)
  F, M = u[1:2], u[3]
  f[4:5] = 1/robot.mass*F
  f[6] = robot.Jinv*M
end

function A_dyn(x::Vector, robot::Robot, model::AstrobeeSE2)
  kron([0 1; 0 0], eye(3))
end

function B_dyn(x::Vector, robot::Robot, model::AstrobeeSE2)
  B = zeros(6,3)
  B[1:6,1:3] = kron([0;1], eye(3))
  B[4:5,:] = 1/robot.mass * B[4:5,:]
  B[6,:] = robot.Jinv * B[6,:]
  return B
end

# Generate full discrete update version of dynamics matrices for a time step
# TODO(ambyld): Rename these? Initialize these once?
function A_dyn_discrete(x, dt, robot::Robot, model::AstrobeeSE2)
  kron([1 dt; 0 1], eye(3))
end

function B_dyn_discrete(x, dt, robot::Robot, model::AstrobeeSE2)
  Jzz_inv = robot.Jinv
  B = [0.5*dt^2*eye(3); dt*eye(3)]
  B[1:2,:] *= 1/robot.mass  # translational force -> acceleration
  B[4:5,:] *= 1/robot.mass
  B[3,:] *= Jzz_inv      # rotational moment -> acceleration
  B[6,:] *= Jzz_inv 
  return B
end

## Convex state inequality constraints
function csi_translational_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  return norm(X[4:5,k]) - robot.hard_limit_vel
end

function csi_angular_velocity_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  return norm(X[6,k]) - robot.hard_limit_omega
end

## Convex control inequality constraints
function cci_translational_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  return 1/robot.mass*norm(U[1:2,k]) - robot.hard_limit_accel
end

function cci_angular_accel_bound(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  Jzz_inv = robot.Jinv
  return abs(Jzz_inv*U[3,k]) - robot.hard_limit_alpha
end

## Nonconvex state inequality constraints
function ncsi_obstacle_avoidance_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)

  rb_idx, env_idx = 1, i
  env_ = WS.btenvironment_keepout

  clearance = model.clearance 

  r = [X[1:2,k];0]
  dist,xbody,xobs = BulletCollision.distance(env_,rb_idx,r,env_idx)

  # See Eq. 12b in "Convex optimization for proximity maneuvering of a spacecraft with a robotic manipulator"
  return clearance - dist
end


## Nonconvex state inequality constraints (convexified)
function ncsi_body_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)

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

function ncsi_arm_obstacle_avoidance_constraints_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)

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
function stri_state_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)

  return norm(X[:,k]-Xp[:,k])^2
end

## Convex state inequality constraints
function ctri_control_trust_region(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int) where {T,E}
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  return norm(U[:,k]-Up[:,k])
end

function get_workspace_location(traj, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, k::Int, i::Int=0) where {T,E}
  return [traj.X[1:2,k]; 0]
end

## Constructing full list of constraint functions
function SCPConstraints(SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
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
    push!(SCPC.convex_state_ineq, (csi_translational_velocity_bound, k, 0))
    push!(SCPC.convex_state_ineq, (csi_angular_velocity_bound, k, 0))
  end

  ## Nonconvex state equality constraints (convexified)
  nothing

  ## Nonconvex state inequality constraints (convexified)
  env_ = WS.btenvironment_keepout
  for k = 1:N, i = 1:length(env_.convex_env_components)
    push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_body_obstacle_avoidance_constraints_convexified, k, i))
    # push!(SCPC.nonconvex_state_convexified_ineq, (ncsi_arm_obstacle_avoidance_constraints_convexified, k, i))
  end

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

# TODO: Generalize this? Need to make A always a vector
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
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

function trust_region_ratio_trajopt(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
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

function trust_region_ratio_mao(traj, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  num,den = 0, 0

  cost_true_prev = cost_true(traj_prev, traj_prev, SCPP)
  cost_true_new = cost_true(traj, traj, SCPP)
  for k in 1:N-1
    cost_true_prev += norm(dynamics_constraints(traj_prev, traj_prev, SCPP, k, 0), Inf)
    cost_true_new += norm(dynamics_constraints(traj, traj, SCPP, k, 0), Inf)
  end
  cost_linearized_new = cost_true(traj, traj, SCPP)
  return (cost_true_prev-cost_true_new)/(cost_true_prev-cost_linearized_new)
end

function get_dual_cvx(prob::Convex.Problem, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}, solver) where {T,E}
	if solver == "Mosek"
		return -MathProgBase.getdual(prob.model)[1:SCPP.PD.model.x_dim]
	else
		return []
	end
end

function get_dual_jump(SCPS::SCPSolution, SCPP::SCPProblem{Astrobee2D{T}, AstrobeeSE2, E}) where {T,E}
	@show -MathProgBase.getconstrduals(SCPS.solver_model.internalModel)[1:SCPP.PD.model.x_dim]
end