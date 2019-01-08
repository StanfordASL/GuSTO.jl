
include("dynamics/blankdynamics.jl")
include("dynamics/freeflyer_se2.jl")
include("dynamics/astrobee_se2.jl")
include("dynamics/astrobee_se3.jl")
include("dynamics/astrobee_se3_manifold.jl")
include("dynamics/airplane.jl")
include("dynamics/dubins_car.jl")
include("dynamics/doubleint_r2.jl")
include("dynamics/panda_kinematics.jl")

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
