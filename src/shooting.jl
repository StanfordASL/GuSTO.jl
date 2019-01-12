export solve!

# TODO: add ability to enforce a goal on part of the states
function solve!(SS::ShootingSolution, SP::ShootingProblem)
	model, x_init, x_goal = SP.PD.model, SP.PD.x_init, SP.PD.x_goal

	# Set up shooting function
	shooting_eval! = (F, p0) -> parameterized_shooting_eval!(F, p0, SP)

	# Run Newton method
	sol_newton = nlsolve(shooting_eval!, SP.p0, iterations = 20) # TODO(ambyld): Make number of iterations a parameter

	if sol_newton.f_converged
		# Recover trajectory
		N, tf, x_dim = SP.N, SP.tf, model.x_dim
		x0 = [x_init; sol_newton.zero]
		tspan = (0., tf)
		dt = tf/(N-1)
		prob = ODEProblem(model_ode!, x0, tspan, SP)
		sol_ode = DifferentialEquations.solve(prob, saveat=dt)
		xp = hcat(sol_ode.u...)
		X, P = xp[1:x_dim,:], xp[x_dim+1:end,:]
		U = get_control(X, P, SP)
		new_traj = Trajectory(X, U, tf, dt)

		# TODO: Add metrics, check for convergence over multiple successes
		push!(SS.prob_status, :Optimal)
	else
		# TODO: Add metrics
		push!(SS.prob_status, :Diverged)
	end
end

function parameterized_shooting_eval!(F, p0, SP::ShootingProblem)
	model, x_init, x_goal = SP.PD.model, SP.PD.x_init, SP.PD.x_goal
	x_dim = model.x_dim
	tf = SP.tf

	x0 = [x_init; p0]
	tspan = (0., tf)
	prob = ODEProblem(model_ode!, x0, tspan, SP)
	sol = DifferentialEquations.solve(prob)

	for i = 1:x_dim
		F[i] = x_goal[i] - sol.u[end][i]
	end
end