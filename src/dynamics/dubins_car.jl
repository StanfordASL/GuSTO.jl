export DubinsCar
export init_traj_nothing, init_traj_straightline

mutable struct DubinsCar <: DynamicsModel
	x_dim	# x, y, theta
	u_dim	# turning rate control

	v 		# speed
	k 		# curvature parameter

	x_min
	x_max
	u_min
	u_max
	clearance

	f::Vector
	A::Vector
	B
end

function DubinsCar()
	x_dim = 3
	u_dim = 1
	v = 2.
	k = 1.
	x_max = [100., 100., 2*π]
	x_min = -x_max
	u_max = 10.
	u_min = -u_max
	clearance = 0.01
	DubinsCar(x_dim, u_dim, v, k, x_min, x_max, u_min, u_max, clearance, [], [], [])
end

function SCPParam(model::DubinsCar, fixed_final_time::Bool)
  convergence_threshold = 0.1

  SCPParam(fixed_final_time, convergence_threshold)
end

function SCPParam_GuSTO(model::DubinsCar)
  Δ0 = 10000.
  ω0 = 1.
  ω_max = 1.0e10
  ε = 1.0e-6
  ρ0 = 0.4
  ρ1 = 1.5 
  β_succ = 2.
  β_fail = 0.5
  γ_fail = 5.

  SCPParam_GuSTO(Δ0, ω0, ω_max, ε, ρ0, ρ1, β_succ, β_fail, γ_fail)
end

###########
# Convex.jl 
###########
function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
	U = traj.U
  N = SCPP.N
  dtp = traj_prev.dt 	# TODO(ambyld): Change to current trajectory dt?

  return dtp*sum(U[j]^2 for j = 1:N-1)
end

function cost_true_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
	cost_true(traj, traj_prev, SCPP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_nothing(TOP::TrajectoryOptimizationProblem{Car, DubinsCar, E}) where E
	model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
	x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

	X = repmat(0.5(x_init + x_goal),1,N)
	U = zeros(u_dim,N-1)
	Trajectory(X, U, tf_guess)
end

function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Car, DubinsCar, E}) where E
	model, x_init, x_goal = TOP.PD.model, TOP.PD.x_init, TOP.PD.x_goal
	x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

	X = hcat(range(x_init, stop=x_goal, length=N)...)
	U = zeros(u_dim,N-1)
	Trajectory(X, U, tf_guess)
end

####################
# Constraint-adding 
####################
function initialize_model_params!(SCPP::SCPProblem{Car, DubinsCar, E}, traj_prev::Trajectory) where E
	N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
	x_dim, u_dim = model.x_dim, model.u_dim
	Xp, Up = traj_prev.X, traj_prev.U

	model.f, model.A, model.B = [], [], B_dyn(Xp[:,1],robot,model)
	for k = 1:N-1
		push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
		push!(model.A, A_dyn(Xp[:,k],robot,model))
	end
end

function update_model_params!(SCPP::SCPProblem{Car, DubinsCar, E}, traj_prev::Trajectory) where E
	N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
	X, U, f, A = traj_prev.X, traj_prev.U, model.f, model.A

	for k = 1:N-1
		update_f!(f[k], X[:,k], U[:,k], robot, model)
		update_A!(A[k], X[:,k], robot, model)
	end
end

macro constraint_abbrev_dubinscar(traj, traj_prev, SCPP)
	quote
		X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
	end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}, k::Int) where E
	# Where i is the state index, and k is the timestep index
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_dubinscar(traj, traj_prev, SCPP)

	fp, Ap, Bp = get_f(k, model), get_A(k, model), get_B(k, model)
  # Tf.*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp*(U[k]-Up[k])) - (X[:,k+1]-X[:,k])/dh

  if k == N-1
    return Tf.*fp + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp.*(U[:,k]-Up[:,k])) - (X[:,k+1]-X[:,k])/dh
  else
    return 0.5*(Tf.*(fp + get_f(k+1, model))  + Tfp*(Ap*(X[:,k]-Xp[:,k]) + Bp.*(U[:,k]-Up[:,k])) +
      Tfp*(Ap*(X[:,k+1]-Xp[:,k+1]) + Bp.*(U[:,k+1]-Up[:,k+1]))) - (X[:,k+1]-X[:,k])/dh
  end
end

# Get current dynamics structures for a particular time step
get_f(k::Int, model::DubinsCar) = model.f[k]
get_A(k::Int, model::DubinsCar) = model.A[k]
get_B(k::Int, model::DubinsCar) = model.B

# Generate full dynamics structures for one time step
function f_dyn(x::Vector, u::Vector, robot::Robot, model::DubinsCar)
	x_dim = model.x_dim
	f = zeros(x_dim)
  update_f!(f, x, u, robot, model)
  return f
end

function update_f!(f, x::Vector, u::Vector, robot::Robot, model::DubinsCar)
	f[1] = model.v*cos(x[3])
	f[2] = model.v*sin(x[3])
	f[3] = model.k*u[1]
end

function A_dyn(x::Vector, robot::Robot, model::DubinsCar)
	x_dim = model.x_dim
	A = zeros(x_dim, x_dim)
  update_A!(A, x, robot, model)
  return A
end

function update_A!(A, x::Vector, robot::Robot, model::DubinsCar)
	A[1,3] = -model.v*sin(x[3])
	A[2,3] = model.v*cos(x[3])
end

function B_dyn(x::Vector, robot::Robot, model::DubinsCar)
  [0; 0; model.k]
end

function add_constraint_category!(SCPCCat, func, dimtype, ind_time, ind_other...)
	SCPCCat[Symbol(func)] = ConstraintCategory(func, dimtype, ind_time, ind_other)
end

## Constructing full list of constraint functions
function SCPConstraints(SCPP::SCPProblem{Car, DubinsCar, E}) where E
	model = SCPP.PD.model
	x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

	SCPC = SCPConstraints()

	## Convex state equality constraints
	# Init and goal (add init first for convenience of getting dual)
	add_constraint_category!(SCPC.convex_state_eq, cse_init_constraints, :scalar, 0, 1:x_dim)
	add_constraint_category!(SCPC.convex_state_eq, cse_goal_constraints, :scalar, 0, 1:x_dim)

  ## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 1:N-1)

	## Convex state inequality constraints
	# State bounds
	add_constraint_category!(SCPC.convex_state_ineq, csi_max_bound_constraints, :scalar, 1:N, 1:x_dim)
	add_constraint_category!(SCPC.convex_state_ineq, csi_min_bound_constraints, :scalar, 1:N, 1:x_dim)

	## Nonconvex state equality constraints
	## Nonconvex state inequality constraints
	## Nonconvex state equality constraints (convexified)
	## Nonconvex state inequality constraints (convexified)
	## Convex control equality constraints
	nothing

	## Convex control inequality constraints
	# Control bounds
	add_constraint_category!(SCPC.convex_control_ineq, cci_max_bound_constraints, :scalar, 1:N-1, 1)
	add_constraint_category!(SCPC.convex_control_ineq, cci_min_bound_constraints, :scalar, 1:N-1, 1)

	SCPC
end

# TODO(ambyld): Replace f_dyn
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
  # Where i is the state index, and k is the timestep index
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev_astrobeeSE2(traj, traj_prev, SCPP)
  fp, Ap = model.f, model.A
  num, den = 0, 0

  for k in 1:N-1
    linearized = fp[k] + Ap[k]*(X[:,k]-Xp[:,k])
    num += norm(f_dyn(X[:,k],U[:,k],robot,model) - linearized)
    den += norm(linearized)
  end

  return num/den
end

##################
# Shooting Method
##################
function get_dual_cvx(prob::Convex.Problem, SCPP::SCPProblem{Car, DubinsCar, E}, solver) where E
	if solver == "Mosek"
		return -MathProgBase.getdual(prob.model)[1:SCPP.PD.model.x_dim]
	else
		return []
	end
end

function get_dual_jump(SCPC::SCPConstraints, SCPP::SCPProblem{Car, DubinsCar, E}) where E
	JuMP.dual.([SCPC.convex_state_eq[:cse_init_constraints].con_reference[0,(i,)] for i = SCPC.convex_state_eq[:cse_init_constraints].ind_other[1]])
end

function model_ode!(xdot, x, SP::ShootingProblem{Car, DubinsCar, E}, t) where E
	x, y, th,	px, py, pth = x
	robot, model = SP.PD.robot, SP.PD.model

	# State variables
	xdot[1] = model.v*cos(th)
	xdot[2] = model.v*sin(th)
	xdot[3] = model.k^2*pth/2

	# Dual variables
	xdot[4] = 0
	xdot[5] = 0
	xdot[6] = px*model.v*sin(th) - py*model.v*cos(th)
end

#######
# JuMP
#######
function add_objective!(solver_model::Model, SCPV::SCPVariables, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
  robot, model = SCPP.PD.robot, SCPP.PD.model
  x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

  U = SCPV.U
  N, dt = SCPP.N, SCPP.tf_guess/SCPP.N

  @NLobjective(solver_model, Min, sum(dt*U[i,k]^2 for i = 1:u_dim, k = 1:N-1))
end



# function add_constraints!(solver_model::Model, SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Astrobee3D{T}, AstrobeeSE3Manifold, E}) where {T,E}
#   X, U = SCPV.X, SCPV.U
#   Xp, Up, Tfp, dtp = traj_prev.X, traj_prev.U, traj_prev.Tf, traj_prev.dt
#   Tf = copy(Tfp)
#   robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
#   x_dim, u_dim, N, dh = model.x_dim, model.u_dim, SCPP.N, SCPP.dh
#   WS = SCPP.WS

# function add_constraints!(solver_model::Model, SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
# 	X, U = SCPV.X, SCPV.U
# 	Xp, Up, dtp = traj_prev.X, traj_prev.U, traj_prev.dt
# 	robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
# 	x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

# 	# Init (add these first for shooting method!)
# 	for i = 1:x_dim
# 		@constraint(solver_model, X[i,1] == x_init[i])
# 	end

# 	# Goal
# 	for i = 1:x_dim
# 		@constraint(solver_model, X[i,N] == x_goal[i])
# 	end

# 	# Dynamics
# 	for k = 1:N-1
# 		@constraint(solver_model, X[1,k+1] == X[1,k] + dtp*model.v*(cos(Xp[3,k]) - sin(Xp[3,k])*(X[3,k] - Xp[3,k])))
# 		@constraint(solver_model, X[2,k+1] == X[2,k] + dtp*model.v*(sin(Xp[3,k]) + cos(Xp[3,k])*(X[3,k] - Xp[3,k])))
# 		@constraint(solver_model, X[3,k+1] == X[3,k] + dtp*model.k*U[k])
# 	end
# end

# function add_nonlinear_constraints!(solver_model::Model, SCPV::SCPVariables, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
# 	X, U = SCPV.X, SCPV.U
# 	Xp, Up, dtp = traj_prev.X, traj_prev.U, traj_prev.dt
# 	robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
# 	x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

# 	# Dynamics
# 	for k = 1:N-1
# 		@NLconstraint(solver_model, X[1,k+1] == X[1,k] + dt*model.v*cos(X[3,k]))
# 		@NLconstraint(solver_model, X[2,k+1] == X[2,k] + dt*model.v*sin(X[3,k]))
# 		@constraint(solver_model, X[3,k+1] == X[3,k] + dt*model.k*U[k])
# 	end

# 	# Init and goal
# 	for i = 1:x_dim
# 		@constraint(solver_model, X[i,1] == x_init[i])
# 		@constraint(solver_model, X[i,N] == x_goal[i])
# 	end
# end
