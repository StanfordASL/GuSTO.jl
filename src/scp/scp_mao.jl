
export solve_mao_cvx!

mutable struct SCPParam_Mao <: SCPParamSpecial
  rho::Vector         # trust-region acceptance thresholds
  Delta_u0            # initial trust region size for control
  lambda              # penalty weight
  alpha               # trust region shrinkage & growth factor 

  r_vec::Vector           # trust region ratios
  Delta_u_vec::Vector     # trust region size for control
end


function SCPParam_Mao(rho, Delta_u0, lambda, alpha)
	SCPParam_Mao(rho, Delta_u0, lambda, alpha, [0.], [Delta_u0])
end

function solve_mao_cvx!(SCPS::SCPSolution, SCPP::SCPProblem, solver="Mosek", max_iter=50, force=false; kwarg...)
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
	param = SCPP.param

	param.alg = SCPParam_Mao(SCPP.PD.model)
  alpha, Delta_u0 = param.alg.alpha, param.alg.Delta_u0
  r_vec, Delta_u_vec = param.alg.r_vec, param.alg.Delta_u_vec
  rho0, rho1, rho2 = param.alg.rho

	iter_cap = SCPS.iterations + max_iter
	SCPV = SCPVariables{Convex.Variable,Convex.Variable}(SCPP)
	SCPC = SCPConstraints(SCPP)

	update_model_params!(SCPP, SCPS.traj)
	push!(SCPS.J_true, cost_true(SCPS.traj, SCPS.traj, SCPP))
	param.obstacle_toggle_distance = model.clearance + 1. # TODO(ambyld): Generalize this, workspace-dependent

	first_time = true
  prob = minimize(0.)
	while (SCPS.iterations <= iter_cap)
		tic()

		# Set up, solve problem
		prob.objective = cost_full_convexified_mao(SCPV, SCPS.traj, SCPC, SCPP)
		prob.constraints = add_constraints_mao_cvx(SCPV, SCPS.traj, SCPC, SCPP)
    Convex.solve!(prob,warmstart=!first_time)
    first_time = false

		push!(SCPS.prob_status, prob.status)
    if prob.status != :Optimal
      warn("SCP-Mao failed find optimal solution")
		  push!(SCPS.iter_elapsed_times,toq()) 
      return
    end
  
		# Recover solution
    new_traj = Trajectory(SCPV.X.value, SCPV.U.value, SCPV.Tf.value[1])
		push!(SCPS.convergence_measure, convergence_metric(new_traj, SCPS.traj, SCPP))
    push!(SCPS.J_full, prob.optval)
		SCPS.dual = get_dual_cvx(prob, SCPP, solver)

    # Check trust regions
		push!(r_vec, trust_region_ratio_mao(new_traj, SCPS.traj, SCPP))
    if r_vec[end] < rho0
			push!(SCPS.accept_solution, false)
      push!(SCPS.J_true, SCPS.J_true[end])
    else
			push!(SCPS.accept_solution, true)
			push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
			copy!(SCPS.traj, new_traj)
    end

    # Adjust trust region radius
    if r_vec[end] < rho0 || r_vec[end] < rho1
      # shrink trust region
      push!(Delta_u_vec, 1/alpha*Delta_u_vec[end])
    elseif r_vec[end] < rho2
      # maintain trust region
      push!(Delta_u_vec, Delta_u_vec[end])
    else
      # expand trust region
      push!(Delta_u_vec, alpha*Delta_u_vec[end])
    end

    iter_elapsed_time = toq()
    push!(SCPS.iter_elapsed_times, iter_elapsed_time)
    SCPS.total_time += iter_elapsed_time
    SCPS.iterations += 1

    !SCPS.accept_solution[end] && continue

		if SCPS.convergence_measure[end] <= param.convergence_threshold
			SCPS.converged = true
			force ? continue : break
		end
	end
end

function add_constraints_mao_cvx(SCPV::SCPVariables, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	constraints = Convex.Constraint[]
	Delta_u = SCPP.param.alg.Delta_u_vec[end]

	update_model_params!(SCPP, traj_prev)
  
  for (f, k, i) in SCPC.convex_state_eq
		constraints += f(SCPV, traj_prev, SCPP, k, i) == 0.
	end

  for (f, k, i) in (SCPC.convex_state_ineq..., SCPC.convex_control_ineq...)
		constraints += f(SCPV, traj_prev, SCPP, k, i) <= 0.
	end

  for (f, k, i) in SCPC.nonconvex_state_convexified_ineq
    # hack for enforcing SDF constraint when it passes out
    # 0 (i.e. too far from obstacle) rather than Convex expr
    new_constraint = f(SCPV, traj_prev, SCPP, k, i)
    typeof(new_constraint <= 0) == Bool && continue
    constraints += new_constraint <= 0.
	end

	for (f, k, i) in SCPC.control_trust_region_ineq
		constraints += f(SCPV, traj_prev, SCPP, k, i) - Delta_u <= 0.
	end

	return constraints
end

function cost_nonconvex_penalty_convexified_mao(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	lambda = SCPP.param.alg.lambda

	J = 0
  for (f, k, i) in SCPC.dynamics
		J += lambda*norm(f(traj, traj_prev, SCPP, k, i), Inf)
	end

	return J
end

function cost_full_convexified_mao(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_true_convexified(traj, traj_prev, SCPP) + cost_nonconvex_penalty_convexified_mao(traj, traj_prev, SCPC, SCPP)
end
