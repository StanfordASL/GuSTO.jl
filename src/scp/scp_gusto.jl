
export solve_gusto_cvx!, solve_gusto_jump!#, solve_gusto_jump_nonlinear!

mutable struct SCPParam_GuSTO <: SCPParamSpecial
	Delta0 								# trust region size 
	omega0         				# exact penalty violation
	omegamax
	epsilon
	rho0
	rho1
	# alpha_succ
	# alpha_fail
  beta_succ		   				# trust-region growth factors
  beta_fail
  gamma_fail

  Delta_vec::Vector
  omega_vec::Vector
  rho_vec::Vector           # trust region ratios
  trust_region_satisfied_vec::Vector  # Bool flags
  convex_ineq_satisfied_vec::Vector  # Bool flags
end

function SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail)
	SCPParam_GuSTO(Delta0, omega0, omegamax, epsilon, rho0, rho1, beta_succ, beta_fail, gamma_fail, [Delta0], [omega0], [0.], [false], [false])
end

######
# CVX 
######

function solve_gusto_cvx!(SCPS::SCPSolution, SCPP::SCPProblem, solver="Mosek", max_iter=50, force=false; kwarg...)
	# Solves a sequential convex programming problem
	# Inputs:
	#		SCPS - SCP solution structure, including an initial trajectory
	#		SCPP - SCP problem definition structure
	#   max_iter - Set to enforce a maximum number of iterations (0 = unlimited)
	#   force - Set true to force further refinement iterations despite convergence

	if solver == "Mosek"
		set_default_solver(MosekSolver(; kwarg...))
	elseif solver == "Gurobi"
		set_default_solver(GurobiSolver(; kwarg...))
  else
		set_default_solver(SCSSolver(; kwarg...))
	end
	
	N = SCPP.N
	param, model = SCPP.param, SCPP.PD.model

	param.alg = SCPParam_GuSTO(model)
	Delta0, omega0, omegamax, rho0, rho1 = param.alg.Delta0, param.alg.omega0, param.alg.omegamax, param.alg.rho0, param.alg.rho1
	beta_succ, beta_fail, gamma_fail = param.alg.beta_succ, param.alg.beta_fail, param.alg.gamma_fail
	Delta_vec, omega_vec, rho_vec = param.alg.Delta_vec, param.alg.omega_vec, param.alg.rho_vec
  trust_region_satisfied_vec, convex_ineq_satisfied_vec = param.alg.trust_region_satisfied_vec, param.alg.convex_ineq_satisfied_vec


	iter_cap = SCPS.iterations + max_iter
	SCPV = SCPVariables{Convex.Variable,Convex.Variable}(SCPP)
	SCPC = SCPConstraints(SCPP)

	initialize_model_params!(SCPP, SCPS.traj)
  push!(SCPS.J_true, cost_true(SCPS.traj, SCPS.traj, SCPP))
	push!(rho_vec, trust_region_ratio_gusto(SCPS.traj, SCPS.traj, SCPP))
	param.obstacle_toggle_distance = Delta_vec[end]/8 + model.clearance # TODO: Generalize clearance

  first_time = true
  prob = minimize(0.)
	while (SCPS.iterations < iter_cap)
		tic()

		# Set up, solve problem
		update_model_params!(SCPP, SCPS.traj)
		prob.objective = cost_full_convexified_gusto(SCPV, SCPS.traj, SCPC, SCPP)
		prob.constraints = add_constraints_gusto_cvx(SCPV, SCPS.traj, SCPC, SCPP)
		Convex.solve!(prob, warmstart=!first_time)
		first_time = false

		push!(SCPS.prob_status, prob.status)
    if prob.status != :Optimal
      warn("GuSTO SCP iteration failed to find an optimal solution")
      push!(SCPS.iter_elapsed_times,toq()) 
      return
    end

		# Recover solution
		new_traj = Trajectory(SCPV.X.value, SCPV.U.value, SCPV.Tf.value[1])
		push!(SCPS.convergence_measure, convergence_metric(new_traj, SCPS.traj, SCPP))
		push!(SCPS.J_full, prob.optval)
		SCPS.dual = get_dual_cvx(prob, SCPP, solver)
		
		# Check trust regions
		push!(trust_region_satisfied_vec, trust_region_satisfied_gusto(new_traj, SCPS.traj, SCPP))
		push!(convex_ineq_satisfied_vec, convex_ineq_satisfied_gusto(new_traj, SCPS.traj, SCPC, SCPP))

    if trust_region_satisfied_vec[end]
    	push!(rho_vec, trust_region_ratio_gusto(new_traj, SCPS.traj, SCPP))			
			if rho_vec[end] > rho1
				push!(SCPS.accept_solution, false)
				push!(Delta_vec, beta_fail*Delta_vec[end])
				push!(omega_vec, omega_vec[end])
			else
				push!(SCPS.accept_solution, true)
				rho_vec[end] < rho0 ? push!(Delta_vec, min(beta_succ*Delta_vec[end], Delta0)) : push!(Delta_vec, Delta_vec[end])
				!convex_ineq_satisfied_vec[end] ? push!(omega_vec, gamma_fail*omega_vec[end]) : push!(omega_vec, omega0)
			end
		else
			push!(SCPS.accept_solution, false)
			push!(Delta_vec, Delta_vec[end])
			push!(omega_vec, gamma_fail*omega_vec[end])
		end

		if SCPS.accept_solution[end]
			push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
			copy!(SCPS.traj, new_traj)	# TODO(ambyld): Maybe deepcopy not needed?
		else
			push!(SCPS.J_true, SCPS.J_true[end])
		end
		param.obstacle_toggle_distance = Delta_vec[end]/8 + model.clearance

		iter_elapsed_time = toq()
		push!(SCPS.iter_elapsed_times, iter_elapsed_time)
		SCPS.total_time += iter_elapsed_time
		SCPS.iterations += 1

		if omega_vec[end] > omegamax
			warn("GuSTO SCP omegamax exceeded")
			break
		end
		!SCPS.accept_solution[end] ? continue : nothing

		if SCPS.convergence_measure[end] <= param.convergence_threshold
			SCPS.converged = true
			convex_ineq_satisfied_vec[end] && (SCPS.successful = true)
			force ? continue : break
		end
	end

end

# Define this for a dynamics model if model parameters need to be updated
# at each SCP iteration
update_model_params!(SCPP::SCPProblem, traj_prev::Trajectory) = nothing

function add_constraints_gusto_cvx(SCPV::SCPVariables, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	constraints = Convex.Constraint[]

	update_model_params!(SCPP, traj_prev)

	for (f, k, j, i) in (SCPC.convex_state_eq..., SCPC.dynamics...)
		constraints += f(SCPV, traj_prev, SCPP, k, j, i) == 0.
	end
	
  for (f, k, j, i) in SCPC.convex_state_goal_ineq
		constraints += f(SCPV, traj_prev, SCPP, k, j, i) <= 0.
	end

  for (f, k, j, i) in SCPC.convex_control_ineq
		constraints += f(SCPV, traj_prev, SCPP, k, j, i) <= 0.
	end

	if !SCPP.param.fixed_final_time
		constraints += SCPV.Tf > 0.1
	end

	return constraints
end

function cost_convex_penalty_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	omega, Delta = SCPP.param.alg.omega_vec[end], SCPP.param.alg.Delta_vec[end]
	J = 0.
	for (f, k, j, i) in SCPC.convex_state_ineq
		J += omega*max(f(traj, traj_prev, SCPP, k, j, i), 0.)
	end
	for (f, k, i) in SCPC.state_trust_region_ineq
		J += omega*max(f(traj, traj_prev, SCPP, k, i) - Delta, 0.)
	end
	return J
end

function cost_nonconvex_penalty_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	omega = SCPP.param.alg.omega_vec[end]
	J = 0.
	for (f, k, j, i) in SCPC.nonconvex_state_ineq
		J += omega*max(f(traj, traj_prev, SCPP, k, j, i), 0.)
	end
	for (f, k, j, i) in SCPC.nonconvex_state_eq
		J += omega*norm(f(traj, traj_prev, SCPP, k, j, i), 1)
	end
	return J
end

function cost_nonconvex_penalty_convexified_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	omega = SCPP.param.alg.omega_vec[end]
	J = 0.
	for (f, k, j, i) in SCPC.nonconvex_state_convexified_ineq
		J += omega*max(f(traj, traj_prev, SCPP, k, j, i), 0.)
	end
	for (f, k, j, i) in SCPC.nonconvex_state_convexified_eq
		J += omega*norm(f(traj, traj_prev, SCPP, k, j, i), 1)
	end
	return J
end

function cost_penalty_full_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_convex_penalty_gusto(traj, traj_prev, SCPC, SCPP) + cost_nonconvex_penalty_gusto(traj, traj_prev, SCPC, SCPP)
end

function cost_penalty_full_convexified_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_convex_penalty_gusto(traj, traj_prev, SCPC, SCPP) + cost_nonconvex_penalty_convexified_gusto(traj, traj_prev, SCPC, SCPP)
end

function cost_full_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
  cost_true(traj, traj_prev, SCPP) + cost_penalty_full_gusto(traj, traj_prev, SCPC, SCPP)
end

function cost_full_convexified_gusto(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_true_convexified(traj, traj_prev, SCPP) + cost_penalty_full_convexified_gusto(traj, traj_prev, SCPC, SCPP)
end

function trust_region_satisfied_gusto(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  # checks for satisfaction of state trust region constraint
  Delta = SCPP.param.alg.Delta_vec[end]
  max_val = -Inf

  for k in 1:SCPP.N
    # norm(traj.X[:,k]-traj_prev.X[:,k])^2 - Delta >= 0 && return false
    val = norm(traj.X[:,k]-traj_prev.X[:,k])^2
    max_val = val > max_val ? val : max_val
  end
  return max_val-Delta <= 0
end

function convex_ineq_satisfied_gusto(traj::Trajectory, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
  # checks for satisfaction of convex state inequalities and nonconvex->convexified state inequalities
	for (f, k, j, i) in (SCPC.convex_state_ineq..., SCPC.nonconvex_state_convexified_ineq...)
		if f(traj, traj_prev, SCPP, k, j, i) > SCPP.param.alg.epsilon
			return false
		end
	end
  return true
end


########
# IPOPT
########
function solve_gusto_jump!(SCPS::SCPSolution, SCPP::SCPProblem, solver="IPOPT", max_iter=Inf, force=false; kwarg...)
	# Solves a sequential convex programming problem using the JuMP interface
	# Inputs:
	#		SCPS - SCP solution structure, including an initial trajectory
	#		SCPP - SCP problem definition structure
	#   max_iter - Set to enforce a maximum number of iterations (0 = unlimited)
	#   force - Set true to force further refinement iterations despite convergence

	# if solver == "IPOPT"
	# 	chosen_solver = IpoptSolver(; kwarg...)
	# end

	N = SCPP.N
	param = SCPP.param

	iter_cap = SCPS.iterations + max_iter
	param.alg = SCPParam_GuSTO(SCPP.PD.model)
	omega0, omegamax, rho0, rho1 = param.alg.omega0, param.alg.omegamax, param.alg.rho0, param.alg.rho1
	beta_succ, beta_fail, gamma_fail = param.alg.beta_succ, param.alg.beta_fail, param.alg.gamma_fail
	Delta_vec, omega_vec, rho_vec = param.alg.Delta_vec, param.alg.omega_vec, param.alg.rho_vec

	iter_cap = SCPS.iterations + max_iter

	while (SCPS.iterations <= iter_cap)
		tic()
		SCPS.solver_model = Model(solver=IpoptSolver(; kwarg...))
		SCPV = SCPVariables{JuMP.Variable, Array{JuMP.Variable}}()
		add_variables!(SCPS.solver_model, SCPV, SCPP)
		add_objective!(SCPS.solver_model, SCPV, SCPP)
		add_constraints!(SCPS.solver_model, SCPV, SCPS.traj, SCPP)
		
		setvalue.(SCPV.X, SCPS.traj.X)
		setvalue.(SCPV.U, SCPS.traj.U)
		
		JuMP.solve(SCPS.solver_model)

		new_traj = Trajectory(getvalue(SCPV.X), getvalue(SCPV.U), SCPS.traj.Tf)
		push!(SCPS.convergence_measure, convergence_metric(new_traj, SCPS.traj, SCPP))

		SCPS.dual = get_dual_jump(SCPS, SCPP)
		iter_elapsed_time = toq()
		push!(SCPS.iter_elapsed_times, iter_elapsed_time)
		SCPS.total_time += iter_elapsed_time
		SCPS.iterations += 1

		push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
		copy!(SCPS.traj, new_traj)

		if SCPS.convergence_measure[end] <= param.convergence_threshold
			SCPS.converged = true
			force ? continue : break
		end
	end
end
