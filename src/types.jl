export ProblemDefinition, Trajectory
export TrajectoryOptimizationProblem, TrajectoryOptimizationSolution

abstract type Robot end
abstract type DynamicsModel end
abstract type Environment end

mutable struct Workspace
	btenvironment_keepin
	btenvironment_keepout

	# signed_distances::Matrix
end
function Workspace(rb::Robot, env::Environment)
  btenv_keepin = BulletCollision.BulletStaticEnvironment(rb.btCollisionObject, BulletCollision.collision_world(env.worldAABBmin, env.worldAABBmax))
  for zone in env.keepin_zones
    BulletCollision.add_collision_object!(btenv_keepin, BulletCollision.geometry_type_to_BT(zone))
  end

  btenv_keepout = BulletCollision.BulletStaticEnvironment(rb.btCollisionObject, BulletCollision.collision_world(env.worldAABBmin, env.worldAABBmax))
  for zone in (env.keepout_zones..., env.obstacle_set...)
    BulletCollision.add_collision_object!(btenv_keepout, BulletCollision.geometry_type_to_BT(zone))
  end

  Workspace(btenv_keepin, btenv_keepout)
end

mutable struct ProblemDefinition{R<:Robot, D<:DynamicsModel, E<:Environment}
	robot::R
	model::D
	env::E

	x_init
	x_goal
end

mutable struct Trajectory
	X 	# x_dim x N
	U
	Tf
	dt
end

mutable struct TrajectoryOptimizationProblem{R<:Robot, D<:DynamicsModel, E<:Environment}
	PD::ProblemDefinition{R,D,E}
  WS::Workspace
	fixed_final_time::Bool

	N 				# Discretization steps
	tf_guess	# Guess for final time
	dh				# Normalized dt
end
function TrajectoryOptimizationProblem(PD,fixed_final_time,N,tf_guess,dh)
  WS = Workspace(PD.robot,PD.env)
  TrajectoryOptimizationProblem(PD,WS,fixed_final_time,N,tf_guess,dh)
end

abstract type SCPParamSpecial end

mutable struct SCPParam
	fixed_final_time::Bool
	convergence_threshold
	obstacle_toggle_distance

	alg::SCPParamSpecial	# Algorithm-specific parameters

	SCPParam(a::Bool, b) = new(a,b)
end

# Optimization algorithm problem
abstract type OptAlgorithmProblem{R<:Robot, D<:DynamicsModel, E<:Environment} end

mutable struct SCPProblem{R,D,E} <: OptAlgorithmProblem{R,D,E}
	PD::ProblemDefinition{R,D,E}
  WS::Workspace
	param::SCPParam

	N 				# Discretization steps
	tf_guess	# Guess for final time
	dh				# Normalized dt
end

# TODO: Remove Convex.jl
VariableTypes = Union{Convex.Variable, JuMP.VariableRef}

mutable struct SCPVariables{T <: VariableTypes, S<:Union{T,Array{T}}}
	X::S 		# State trajectory
	U::S 		# Control trajectory
	Tf::T 	# Final time

	SCPVariables{T,S}() where {T <: VariableTypes, S<:Union{T,Array{T}}} = new()
	SCPVariables{T,S}(X,U) where {T <: VariableTypes, S<:Union{T,Array{T}}} = new(X,U)
	SCPVariables{T,S}(X,U,Tf) where {T <: VariableTypes, S<:Union{T,Array{T}}} = new(X,U,Tf)
end
SCPVariables(X::S,U::S) where {S<:Union{VariableTypes,Array{VariableTypes}}} = SCPVariables{T,S}(X,U)
SCPVariables(X::S,U::S,Tf::T) where {T <: VariableTypes, S<:Union{T,Array{T}}} = SCPVariables{T,S}(X,U,Tf)
function SCPVariables{T,S}(SCPP::SCPProblem) where {T <: Convex.Variable, S<:Union{T,Array{T}}}
	X = Convex.Variable(SCPP.PD.model.x_dim, SCPP.N)
	U = Convex.Variable(SCPP.PD.model.u_dim, SCPP.N-1)
	Tf = Convex.Variable(1)
	SCPP.param.fixed_final_time ? fix!(Tf, SCPP.tf_guess) : nothing
	SCPVariables{T,S}(X,U,Tf)
end

ProbStatusTypes = Union{Symbol, MOI.TerminationStatusCode}

mutable struct SCPConstraints
	dynamics::Dict
	convex_state_eq::Dict
	convex_state_ineq::Dict
	nonconvex_state_eq::Dict
	nonconvex_state_ineq::Dict
	nonconvex_state_convexified_eq::Dict
	nonconvex_state_convexified_ineq::Dict
	convex_state_bc_eq::Dict 
	convex_state_bc_ineq::Dict    # TODO(acauligi): only convex_state_bc_ineq, nonconvex_state_bc_ineq, nonconvex_state_bc_convexified_ineq
	nonconvex_state_bc_eq::Dict   #   are correctly implemented across algorithms
	nonconvex_state_bc_ineq::Dict
	nonconvex_state_bc_convexified_eq::Dict
	nonconvex_state_bc_convexified_ineq::Dict
	convex_control_eq::Dict
	convex_control_ineq::Dict
	# obstacle_avoidance_ineq::Vector
	# obstacle_avoidance_convexified_ineq::Vector
	state_trust_region_ineq::Dict
	control_trust_region_ineq::Dict
end
SCPConstraints() = SCPConstraints((Dict{Symbol, ConstraintCategory}() for i in 1:17)...)

mutable struct SCPSolution
	traj::Trajectory
	dual::Vector

	J_true::Vector
	J_full::Vector
	prob_status::Vector{ProbStatusTypes}	# Solver status
	accept_solution::Vector 							# Solution accepted?
	convergence_measure::Vector						# Convergence measure of the latest iteration
	successful::Bool 											# Did we find an acceptable solution?
	converged::Bool												# Has the solution met the convergence condition?
	iterations::Int 											# Number of SCP iterations executed
	iter_elapsed_times::Vector
	total_time 	# TODO(ambyld): Move to TOP

	param::SCPParam
	SCPP::SCPProblem
	SCPC::SCPConstraints
	solver_model

	SCPSolution(a,b,c,d,e,f,g,h,i,j,k,l,m,n) = new(a,b,c,d,e,f,g,h,i,j,k,l,m,n)
end

mutable struct ConstraintCategory
	func
	dimtype
	ind_time
	ind_other
	con_reference
	var_reference
	ConstraintCategory(a,b,c,d) = new(a,b,c,d)
end

mutable struct ShootingProblem{R,D,E} <: OptAlgorithmProblem{R,D,E}
	PD::ProblemDefinition{R,D,E}
	WS::Workspace

	p0 	# Initial dual
	N 	# Discretization steps
	tf  # Final time
end

mutable struct ShootingSolution
	traj::Trajectory

	J_true::Vector
	prob_status::Vector
	convergence_measure::Vector
	converged::Bool
	iter_elapsed_times::Vector

	SP::ShootingProblem
end
ShootingSolution(SP, traj_init) = ShootingSolution(traj_init, [], [:(NA)], [NaN], false, [0.], SP)

mutable struct TrajectoryOptimizationSolution
	traj::Trajectory
	SCPS::SCPSolution
	SS::ShootingSolution
	total_time
	TrajectoryOptimizationSolution(TOP::TrajectoryOptimizationProblem) = new(Trajectory(TOP))
end

ShootingProblem(TOP::TrajectoryOptimizationProblem, SCPS::SCPSolution) = ShootingProblem(TOP.PD, TOP.WS, SCPS.dual, TOP.N, SCPS.traj.Tf)

TrajectoryOptimizationProblem(PD, N, tf_guess; fixed_final_time::Bool=false) = TrajectoryOptimizationProblem(PD, fixed_final_time, N, tf_guess, 1/(N-1))

# Initialize a blank trajectory optimization solution
# TrajectoryOptimizationSolution(TOP::TrajectoryOptimizationProblem) = TrajectoryOptimizationSolution(Trajectory(TOP))

SCPSolution(SCPP::SCPProblem, traj_init::Trajectory) = SCPSolution(traj_init, [], [], [], [:(NA)], [true], [0.], false, false, 0, [0.], 0., SCPP.param, SCPP)

Trajectory(X, U, Tf) = Trajectory(X, U, Tf, Tf/(size(X,2)-1))

# Initialize a blank trajectory structure
function Trajectory(TOP::TrajectoryOptimizationProblem)
	x_dim = TOP.PD.model.x_dim
	u_dim = TOP.PD.model.u_dim
	N, tf_guess = TOP.N, TOP.tf_guess
	dt = tf_guess/(N-1)

	Trajectory(zeros(x_dim,N), zeros(u_dim,N-1), tf_guess, dt)
end

function Base.copy!(a::Trajectory, b::Trajectory)
	a.X = deepcopy(b.X)
	a.U = deepcopy(b.U)
	a.Tf = deepcopy(b.Tf)
	a.dt = deepcopy(b.dt)
end

Base.deepcopy(a::Trajectory) = Trajectory(deepcopy(a.X), deepcopy(a.U), deepcopy(a.Tf), deepcopy(a.dt))

function SCPProblem(TOP::TrajectoryOptimizationProblem)
	N = TOP.N
	SCPProblem(TOP.PD, TOP.WS, SCPParam(TOP.PD.model, TOP.fixed_final_time), N, TOP.tf_guess, 1/(N-1))
end

function SCPParam(model::DynamicsModel, fixed_final_time::Bool)
  convergence_threshold = 0.001
  obstacle_toggle_distance = Inf

  SCPParam(fixed_final_time, convergence_threshold, obstacle_toggle_distance)
end
