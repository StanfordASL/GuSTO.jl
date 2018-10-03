export BlankDynamics

mutable struct BlankDynamics <: DynamicsModel
  x_dim
  u_dim
end

BlankDynamics() = BlankDynamics(3,3)

function cost_true(traj, traj_prev::Trajectory, SCPP::SCPProblem{BlankRobot, BlankDynamics})
	0
end

function cost_true_convexified(traj, traj_prev::Trajectory, SCPP::SCPProblem{BlankRobot, BlankDynamics})
	cost_true(traj, traj_prev, SCPP)
end



macro blank_constraint(expr)
	return quote
		function (traj, traj_prev, SCPP)
			X, U, Tf = traj.X, traj.U, traj.Tf
			Xp, Up, Tfp = traj_prev.X, traj_prev.U, traj_prev.Tf
			robot, model, x_init, x_goal = SCPP.PD.robot, SCPP.PD.model, SCPP.PD.x_init, SCPP.PD.x_goal
			x_dim, u_dim, N = model.x_dim, model.u_dim, SCPP.N

			return (X, U, Tf, Tfp, robot, model, x_init, x_goal, x_dim, u_dim) -> $expr
		end
	end
end

function SCPConstraints(traj, traj_prev::Trajectory, SCPP::SCPProblem{BlankRobot, BlankDynamics})

	SCPC = SCPConstraints()

	# Dynamics constraints
	for k = 1:N-1
		push!(SCPC.dynamics, @blank_constraint X[1,k] + Tf*U[k] - X[1,k+1])
		push!(SCPC.dynamics, @blank_constraint X[2,k] + Tf*U[k] - X[2,k+1])
		push!(SCPC.dynamics, @blank_constraint X[3,k] + Tf*U[k] - X[3,k+1])
	end

	# Convex state constraints

	# Nonconvex state constraints

	# Nonconvex state constraints (convexified)

	# Convex control constraints

	return SCPC
end





