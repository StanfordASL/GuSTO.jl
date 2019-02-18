export solve!

# TODO: add ability to enforce a goal on part of the states
function solve!(SS::ShootingSolution, SP::ShootingProblem)
	model, x_init, x_goal, tf = SP.PD.model, SP.PD.x_init, SP.PD.goal_set, SP.tf, SP.x_goal

	# Set up shooting function
	shooting_eval! = (F, p0) -> parameterized_shooting_eval!(F, p0, SP)

	# Run Newton method
	time_start = time_ns()
	sol_newton = nothing
	try
		sol_newton = nlsolve(shooting_eval!, SP.p0, iterations = 100, ftol=1e-3) # TODO(ambyld): Make number of iterations a parameter
	catch
		iter_elapsed_time = (time_ns() - time_start)/10^9
		push!(SS.prob_status, :Diverged)
		push!(SS.J_true, NaN)
		push!(SS.convergence_measure, NaN)
		push!(SS.iter_elapsed_times, iter_elapsed_time)
		return
	end
	iter_elapsed_time = (time_ns() - time_start)/10^9

	if sol_newton.f_converged
		# Recover trajectory
		N, tf, x_dim = SP.N, SP.tf, model.x_dim
		x0 = [x_init; sol_newton.zero]
		tspan = (0., tf)
		dt = tf/(N-1)
		prob = ODEProblem(shooting_ode!, x0, tspan, SP)
		sol_ode = DifferentialEquations.solve(prob, saveat=dt)
		xp = hcat(sol_ode.u...)
		X, P = xp[1:x_dim,:], xp[x_dim+1:end,:]
		U = get_control(X, P, SP)
		new_traj = Trajectory(X, U, tf, dt)

		push!(SS.prob_status, :Optimal)
		push!(SS.J_true, cost_true(new_traj, new_traj, SP))
		push!(SS.convergence_measure, convergence_metric(new_traj, SS.traj, SP))
		push!(SS.iter_elapsed_times, iter_elapsed_time)
		copy!(SS.traj, new_traj)
	else
		push!(SS.prob_status, :Diverged)
		push!(SS.J_true, NaN)
		push!(SS.convergence_measure, NaN)
		push!(SS.iter_elapsed_times, iter_elapsed_time)
	end
end

function parameterized_shooting_eval!(F::Vector{T}, p0::Vector{T}, SP::ShootingProblem) where T <: AbstractFloat
	model, x_init, x_goal, tf, dt = SP.PD.model, SP.PD.x_init, SP.x_goal, SP.tf, SP.dt
	x_dim = model.x_dim

	x0 = [x_init; p0]
	# N = SP.N
	tspan = (0., tf)
	# dt = tf/(500-1)
	prob = ODEProblem(shooting_ode!, x0, tspan, SP)
	sol = DifferentialEquations.solve(prob, dtmin=dt, force_dtmin=true, saveat=dt)

	a = [1.0*ones(6); 10.0*ones(7)]
	for i = 1:x_dim
		F[i] = x_goal[i] - sol.u[end][i]
	end
end