
include("dynamics/blankdynamics.jl")
include("dynamics/freeflyer_se2.jl")
include("dynamics/astrobee_se2.jl")
include("dynamics/astrobee_se3.jl")
include("dynamics/astrobee_se3_manifold.jl")
include("dynamics/airplane.jl")
include("dynamics/dubins_car.jl")
include("dynamics/doubleint_r2.jl")

macro constraint_shortcut(traj, traj_prev, SCPP)
	quote
		X, U, Tf = $(esc(traj)).X, $(esc(traj)).U, $(esc(traj)).Tf
    Xp, Up, Tfp, dtp = $(esc(traj_prev)).X, $(esc(traj_prev)).U, $(esc(traj_prev)).Tf, $(esc(traj_prev)).dt
    robot, model, WS, x_init, goal_set = $(esc(SCPP)).PD.robot, $(esc(SCPP)).PD.model, $(esc(SCPP)).WS, $(esc(SCPP)).PD.x_init, $(esc(SCPP)).PD.goal_set
    x_dim, u_dim, N, dh = model.x_dim, model.u_dim, $(esc(SCPP)).N, $(esc(SCPP)).dh
    X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh
	end
end

## Dynamics constraints
nothing

## Convex state boundary condition equality constraints
function cse_init_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	X[i,1] - x_init[i]
end

function csbce_goal_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, goal::G) where G <: Goal{PointGoal}
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	ind, point = goal.ind_coordinates, goal.params.point
	# @show ind, point, X[ind,N] - point
	X[ind,N] - point
end

## Convex state boundary condition inequality constraints
function csbci_goal_min_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	X[i,1] - x_init[i]
end

function csbci_goal_max_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	X[i,N] - x_goal[i]
end

## Convex state inequality constraints
function csi_max_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	X[i,k] - model.x_max[i]
end

function csi_min_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	model.x_min[i] - X[i,k]
end

## Nonconvex state equality constraints
## Nonconvex state inequality constraints
## Nonconvex state equality constraints (convexified)
## Nonconvex state inequality constraints (convexified)
## Convex control equality constraints

## Convex control inequality constraints
function cci_max_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	U[i,k] - model.u_max[i]
end

function cci_min_bound_constraints(traj, traj_prev::Trajectory, SCPP::SCPProblem, k::Int, i::Int)
	X,U,Tf,Xp,Up,Tfp,dtp,robot,model,WS,x_init,goal_set,x_dim,u_dim,N,dh = @constraint_shortcut(traj, traj_prev, SCPP)
	model.u_min[i] - U[i,k]
end
