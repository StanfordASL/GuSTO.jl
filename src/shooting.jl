export solve!

function solve!(SS::ShootingSolution, SP::ShootingProblem)
	model, x_init, x_goal = SP.PD.model, SP.PD.x_init, SP.PD.x_goal

	# Set up shooting function
	shooting_eval! = (F, p0) -> parameterized_shooting_eval!(F, p0, SP)

	# Run Newton method
	sol = nlsolve(shooting_eval!, SP.p0, iterations = 20) # TODO(ambyld): Make number of iterations a parameter

	SS.converged = sol.f_converged
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
		F[i] = x_goal[i] - sol[end][i]
	end
end