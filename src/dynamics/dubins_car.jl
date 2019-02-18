export DubinsCar
export init_traj_nothing, init_traj_straightline

mutable struct DubinsCar <: DynamicsModel
	x_dim	# x, y, θ
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
  convergence_threshold = 1e-4
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

function cost_true(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Car, DubinsCar, E} where E
	U, N = traj.U, OAP.N
  dtp = traj_prev.dt 	# TODO(ambyld): Change to current trajectory dt?

  Jm = 0

  # Forward Euler
  # Jm += sum(dt*U[k-1]^2 for k = 2:N)

  # Trapezoidal
  Jm += sum(1/2*dtp*(U[k-1]^2 + U[k]^2) for k = 2:N)
end

function cost_true_convexified(traj, traj_prev::Trajectory, OAP::A) where A <: OptAlgorithmProblem{Car, DubinsCar, E} where E
	cost_true(traj, traj_prev, OAP)
end

#############################
# Trajectory Initializations
#############################
function init_traj_nothing(TOP::TrajectoryOptimizationProblem{Car, DubinsCar, E}) where E
	model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
	x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

	goal_final = get_first_goal_at_time(goal_set, tf_guess)
	x_goal = goal_final.params.point

	X = repmat(0.5(x_init + x_goal),1,N)
	U = zeros(u_dim,N)
	Trajectory(X, U, tf_guess)
end

function init_traj_straightline(TOP::TrajectoryOptimizationProblem{Car, DubinsCar, E}) where E
	model, x_init, goal_set = TOP.PD.model, TOP.PD.x_init, TOP.PD.goal_set
	x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, TOP.N, TOP.tf_guess

	goal_final = get_first_goal_at_time(goal_set, tf_guess)
	x_goal = goal_final.params.point

	X = hcat(range(x_init, stop=x_goal, length=N)...)
	U = zeros(u_dim,N)
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
	for k = 1:N
		push!(model.f, f_dyn(Xp[:,k],Up[:,k],robot,model))
		push!(model.A, A_dyn(Xp[:,k],robot,model))
	end
end

function update_model_params!(SCPP::SCPProblem{Car, DubinsCar, E}, traj_prev::Trajectory) where E
	N, robot, model = SCPP.N, SCPP.PD.robot, SCPP.PD.model
	Xp, Up, f, A = traj_prev.X, traj_prev.U, model.f, model.A

	for k = 1:N
		update_f!(f[k], Xp[:,k], Up[:,k], robot, model)
		update_A!(A[k], Xp[:,k], robot, model)
	end
end

macro constraint_abbrev_dubinscar(traj, traj_prev, SCPP)
	quote
		X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, goal_set = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.goal_set
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh
	end
end

## Dynamics constraints
function dynamics_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}, k::Int) where E
	# Where i is the state index, and k is the timestep index
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_abbrev_dubinscar(traj, traj_prev, SCPP)

  Xkp, Ukp, Xk = X[:,k-1], U[k-1], X[:,k]
  fpkp, Apkp, Bpkp, Xpkp, Upkp = get_f(k-1, model), get_A(k-1, model), get_B(k-1, model), Xp[:,k-1], Up[k-1]

  # Just Trapezoidal rule
  Uk = U[k]
  fpk, Apk, Bpk, Xpk, Upk = get_f(k, model), get_A(k, model), get_B(k, model), Xp[:,k], Up[k]
  return (Xkp-Xk) + 1/2*dtp.*(fpkp + Apkp*(Xkp-Xpkp) + Bpkp.*(Ukp-Upkp) +
                              fpk  + Apk *(Xk-Xpk)   + Bpk .*(Uk-Upk))
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

## Constructing full list of constraint functions
function SCPConstraints(SCPP::SCPProblem{Car, DubinsCar, E}) where E
	model, goal_set = SCPP.PD.model, SCPP.PD.goal_set
	x_dim, u_dim, N, tf_guess = model.x_dim, model.u_dim, SCPP.N, SCPP.tf_guess

	SCPC = SCPConstraints()

	## Dynamics constraints
  add_constraint_category!(SCPC.dynamics, dynamics_constraints, :array, 2:N)

	## Convex state equality constraints
	# Init and goal
	add_constraint_category!(SCPC.convex_state_eq, cse_init_constraints, :scalar, 0, 1:x_dim)
	for goal in values(inclusive(goal_set.goals, searchsortedfirst(goal_set.goals, tf_guess), searchsortedlast(goal_set.goals, tf_guess)))
		if typeof(goal.params) == PointGoal
			add_constraint_category!(SCPC.convex_state_boundary_condition_eq, csbce_goal_constraints, goal, :array)
		else
			add_constraint_category!(SCPC.convex_state_boundary_condition_ineq, csbci_goal_constraints, goal, :array)
		end
	end

	## Convex state inequality constraints
	# State bounds
	add_constraint_category!(SCPC.convex_state_ineq, csi_max_bound_constraints, :scalar, 1:N, 1:x_dim)
	add_constraint_category!(SCPC.convex_state_ineq, csi_min_bound_constraints, :scalar, 1:N, 1:x_dim)

	## Nonconvex state equality constraints
	nothing
	## Nonconvex state inequality constraints
	nothing
	## Nonconvex state equality constraints (convexified)
	nothing
	## Nonconvex state inequality constraints (convexified)
	nothing
	## Convex control equality constraints
	nothing

	## Convex control inequality constraints
	# Control bounds
	add_constraint_category!(SCPC.convex_control_ineq, cci_max_bound_constraints, :scalar, 1:N-1, 1)
	add_constraint_category!(SCPC.convex_control_ineq, cci_min_bound_constraints, :scalar, 1:N-1, 1)

	return SCPC
end

# TODO(ambyld): Replace f_dyn
function trust_region_ratio_gusto(traj, traj_prev::Trajectory, SCPP::SCPProblem{Car, DubinsCar, E}) where E
  X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_abbrev_dubinscar(traj, traj_prev, SCPP)
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
	init_constraint_category = SCPC.convex_state_eq[:cse_init_constraints][1]
  -JuMP.dual.([init_constraint_category.con_reference[0,(i,)] for i = init_constraint_category.ind_other[1]])
end

function shooting_ode!(xdot, x, SP::ShootingProblem{Car, DubinsCar, E}, t) where E
	X, p = x[1:3], x[4:6]
	U = get_control(X, p, SP)[1]
	xx, y, θ = X
  px, py, pθ = p
	robot, model = SP.PD.robot, SP.PD.model

	# State variables 
	xdot[1] = model.v*cos(θ) 	 	# xdot
	xdot[2] = model.v*sin(θ)		# ydot
	xdot[3] = model.k*U 				# θdot

	# Dual variables
	xdot[4] = 0
	xdot[5] = 0
	xdot[6] = px*model.v*sin(θ) - py*model.v*cos(θ)
end

function get_control(x, p, SP::ShootingProblem{Car, DubinsCar, E}) where E
	model = SP.PD.model
	U = model.k/2*p[3,:]
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
