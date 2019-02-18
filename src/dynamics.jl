
include("dynamics/blankdynamics.jl")
include("dynamics/freeflyer_se2.jl")
include("dynamics/astrobee_se2.jl")
include("dynamics/astrobee_se3.jl")
include("dynamics/astrobee_se3_manifold.jl")
include("dynamics/airplane.jl")
include("dynamics/dubins_car.jl")
include("dynamics/doubleint_r2.jl")

export get_dt_from_SCPP, get_params_from_SCPP
export f_dt, f_dt_dx, f_dt_du

macro constraint_abbrev(traj, traj_prev, SCPP)
	quote
		X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, x_goal = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.x_goal
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh
	end
end

# macro make_func(name, code)
# 	quote
# 		function $name(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
# 			X,U,Tf,Xp,Up,Tfp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dt = @constraint_abbrev(traj, traj_prev, SCPP)
# 			$code
# 		end
# 	end
# end

## Convex state equality constraints
function cse_init_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	X[i,1] - x_init[i]
end

function cse_goal_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	X[i,N] - x_goal[i]
end

## Convex state inequality constraints
function csi_max_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	X[i,k] - model.x_max[i]
end

function csi_min_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	model.x_min[i] - X[i,k]
end

## Nonconvex state equality constraints
## Nonconvex state inequality constraints
## Nonconvex state equality constraints (convexified)
## Nonconvex state inequality constraints (convexified)
## Convex control equality constraints

## Convex control inequality constraints
function cci_max_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	U[i,k] - model.u_max[i]
end

function cci_min_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,x_goal,x_dim,u_dim,N,dh = @constraint_abbrev(traj, traj_prev, SCPP)
	model.u_min[i] - U[i,k]
end





##########################
# Chance constrained MPC #
##########################
macro get_dt_from_SCPP(SCPP)
	quote
	    (($(esc(SCPP)).dh)*($(esc(SCPP)).tf_guess))
	end
end
macro get_params_from_SCPP(SCPP)
	quote
	    robot, model  = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model
	    x_dim, u_dim  = model.x_dim, model.u_dim
	    N, dh, tf, dt = $(esc(SCPP)).N, $(esc(SCPP)).dh, $(esc(SCPP)).tf_guess, ($(esc(SCPP)).dh)*($(esc(SCPP)).tf_guess)
	    
	    robot, model, x_dim, u_dim, N, dh, tf, dt
	end
end


# Get an approximation for the next state given by the dynamics \dot{x} = f(x,u) 
function f_dt(x::Vector{T}, u_step::Vector{T}, SCPP::SCPProblem) where T<:AbstractFloat
    return (x + f_dyn(x, u_step, SCPP.PD.robot, SCPP.PD.model) * @get_dt_from_SCPP(SCPP))
end
function f_dt_dx(x::Vector{T}, u_step::Vector{T}, SCPP::SCPProblem) where T<:AbstractFloat
    # \dot{x} = f(x,u) => x_{k+1} \approx x_{k} + f(x,u)*dt    
	return (eye(model.x_dim) + A_dyn(x, u_step, SCPP.PD.robot, SCPP.PD.model) * @get_dt_from_SCPP(SCPP))
end
function f_dt_du(x::Vector{T}, u_step::Vector{T}, SCPP::SCPProblem) where T<:AbstractFloat
    # \dot{x} = f(x,u) => x_{k+1} \approx x_{k} + f(x,u)*dt    
	return (B_dyn(x, u_step, SCPP.PD.robot, SCPP.PD.model) * @get_dt_from_SCPP(SCPP))
end