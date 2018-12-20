
export solve_trajopt_cvx!

mutable struct SCPParam_TrajOpt <: SCPParamSpecial
  mu0         # initial penalty coefficient
  s0          # initial trust region size
  c           # step acceptance parameter
  τ_plus    # trust region expansion factor
  τ_minus   # trust region shrinkage factor
  k           # penalty scaling factor
  ftol        # convergence threshold for merit
  xtol        # convergence threhsold for x
  ctol        # constraint satisfaction threshold
  max_penalty_iteration
  max_convex_iteration
  max_trust_iteration

  ρ_vec::Vector # trust region ratios
  mu_vec::Vector  # penalty coefficient
  s_vec::Vector   # trust region size
  xtol_vec::Vector
  ftol_vec::Vector
  ctol_vec::Vector
end

function SCPParam_TrajOpt(mu0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol,max_penalty_iteration,max_convex_iteration,max_trust_iteration)
	SCPParam_TrajOpt(mu0, s0, c, τ_plus, τ_minus, k, ftol, xtol, ctol, 
    max_penalty_iteration,max_convex_iteration,max_trust_iteration,[0.],[mu0], [s0], [0.], [0.], [0.])
end

#########
# Gurobi
#########
function solve_trajopt_cvx!(SCPS::SCPSolution, SCPP::SCPProblem, solver="Mosek", max_iter=125, force=false; kwarg...)
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

	param.alg = SCPParam_TrajOpt(SCPP.PD.model)
  mu0,s0 = param.alg.mu0, param.alg.s0
  ρ_vec,mu_vec,s_vec = param.alg.ρ_vec, param.alg.mu_vec, param.alg.s_vec
  τ_plus,τ_minus = param.alg.τ_plus, param.alg.τ_minus
  c = param.alg.c
  ftol,xtol,ctol = param.alg.ftol, param.alg.xtol, param.alg.ctol
  ftol_vec, xtol_vec, ctol_vec = param.alg.ftol_vec, param.alg.xtol_vec, param.alg.ctol_vec

  max_penalty_iteration,max_convex_iteration,max_trust_iteration = 
    param.alg.max_penalty_iteration, param.alg.max_convex_iteration, param.alg.max_trust_iteration

	iter_cap = SCPS.iterations + max_iter
	SCPV = SCPVariables{Convex.Variable,Convex.Variable}(SCPP)
	SCPC = SCPConstraints(SCPP)

  update_model_params!(SCPP, SCPS.traj)
	push!(SCPS.J_true, cost_true(SCPS.traj, SCPS.traj, SCPP))
	param.obstacle_toggle_distance = model.clearance + 1. # TODO(ambyld): Generalize this, workspace-dependent
	
	first_time = true
  prob = minimize(0.)
  constraints_satisfied = false
  xtol_satisfied = false

  old_penalty_traj, old_convex_traj = SCPS.traj, SCPS.traj

  for penalty_iteration in 1:max_penalty_iteration
    constraints_satisfied && break
    copy!(old_penalty_traj,SCPS.traj)

    for convex_iteration in 1:max_convex_iteration
      copy!(old_convex_traj,SCPS.traj)

      constraints_satisfied && break
      if xtol_satisfied
        xtol_satisfied = false
        break
      end
      for trust_iteration in 1:max_trust_iteration
        time_start = time_ns()
        # Set up, solve problem
				prob.objective = cost_full_convexified_trajopt(SCPV, old_convex_traj, SCPC, SCPP)
				prob.constraints = add_constraints_trajopt_cvx(SCPV, old_convex_traj, SCPC, SCPP)
				Convex.solve!(prob, warmstart=!first_time)
				first_time = false

		    push!(SCPS.prob_status, prob.status)
        if prob.status != :Optimal
          warn("TrajOpt failed find optimal solution")
		      push!(SCPS.iter_elapsed_times, (time_ns() - time_start)/10^9) 
          return
        end

        # Recover solution
        new_traj = Trajectory(SCPV.X.value, SCPV.U.value, SCPV.Tf.value[1])
        push!(xtol_vec, evaluate_xtol(new_traj, old_convex_traj, SCPP))
        push!(SCPS.convergence_measure, xtol_vec[end])
        push!(SCPS.J_full, prob.optval)

        SCPS.dual = get_dual_cvx(prob, SCPP, solver)

        # grow or shrink trust region 
        push!(ρ_vec, trust_region_ratio_trajopt(new_traj, old_convex_traj, SCPP))
        if ρ_vec[end] > c
          push!(s_vec, τ_plus*s_vec[end])
        else
          push!(s_vec, τ_minus*s_vec[end])
        end
        
        copy!(SCPS.traj, new_traj)
        iter_elapsed_time = (time_ns() - time_start)/10^9
        push!(SCPS.J_true, cost_true(new_traj, SCPS.traj, SCPP))
        push!(SCPS.iter_elapsed_times, iter_elapsed_time)
        SCPS.total_time += iter_elapsed_time
        SCPS.iterations += 1
        if s_vec[end] < xtol
          xtol_satisfied = true
          break
        end
      end
      
      push!(ftol_vec, evaluate_ftol(SCPS.traj, old_convex_traj, SCPP))
      push!(xtol_vec, evaluate_xtol(SCPS.traj, old_convex_traj, SCPP))
      if ftol_vec[end] < ftol || xtol[end] < xtol
        constraints_satisfied = true
        break
      end
    end

    if evaluate_ctol(SCPS.traj, old_penalty_traj, SCPP, SCPC) < ctol
  	 constraints_satisfied = true
     SCPS.converged = true
     break
  	else 
      # increase penalty coefficient 
      push!(mu_vec, mu_vec[end]*SCPP.param.alg.k)
    end
  end
end

function add_constraints_trajopt_cvx(SCPV::SCPVariables, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
  s = SCPP.param.alg.s_vec[end]
	constraints = Convex.Constraint[]

	update_model_params!(SCPP, traj_prev)

  for (f, k, i) in SCPC.state_trust_region_ineq
		constraints += f(SCPV, traj_prev, SCPP, k, i) - s <= 0
	end
  for (f, k, i) in SCPC.convex_state_eq 
		constraints += f(SCPV, traj_prev, SCPP, k, i) == 0
	end
 
	return constraints
end

function cost_convex_state_eq_penalty_trajopt(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	mu = SCPP.param.alg.mu_vec[end]

	J = 0
	for (f, k, i) in (SCPC.nonconvex_state_convexified_eq..., SCPC.dynamics...)
		J += mu*norm(f(traj, traj_prev, SCPP, k, i),1)
	end
  return J
end

function cost_convex_ineq_penalty_trajopt(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	mu = SCPP.param.alg.mu_vec[end]

	J = 0
	for (f, k, i) in (SCPC.convex_state_ineq..., SCPC.nonconvex_state_convexified_ineq...)
		J += mu*max(f(traj, traj_prev, SCPP, k, i), 0)
  end

  for (f, k, i) in SCPC.convex_control_ineq
		J += mu*max(f(traj, traj_prev, SCPP, k, i), 0)
  end

  return J
end

function cost_penalty_full_convexified_trajopt(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_convex_state_eq_penalty_trajopt(traj, traj_prev, SCPC, SCPP) + cost_convex_ineq_penalty_trajopt(traj, traj_prev, SCPC, SCPP)
end

function cost_full_convexified_trajopt(traj, traj_prev::Trajectory, SCPC::SCPConstraints, SCPP::SCPProblem)
	cost_true_convexified(traj, traj_prev, SCPP) + cost_penalty_full_convexified_trajopt(traj, traj_prev, SCPC, SCPP)
end

function evaluate_xtol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  return convergence_metric(traj,traj_prev,SCPP)
end

function evaluate_ftol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem)
  abs(cost_true(traj, traj, SCPP)-cost_true(traj_prev, traj_prev, SCPP))/abs(cost_true(traj, traj, SCPP))
end

function evaluate_ctol(traj::Trajectory, traj_prev::Trajectory, SCPP::SCPProblem, SCPC::SCPConstraints)
	JNum = 0
  JDen = 0
	constraint_list = (SCPC.convex_state_ineq..., SCPC.nonconvex_state_ineq...,SCPC.convex_state_eq...,SCPC.dynamics...,SCPC.nonconvex_state_eq...)
	if length(constraint_list) > 0
		(f, k, i) = constraint_list[1] 	# Get first function handle
		JtestNum = 0
    JtestDen = 0
		for (fnext, k, i) in constraint_list	# Loop through all constraints
			if fnext == f	# Evaluate and update maximum value while considering same class of constraints
				JtestNum = max(JtestNum, norm(fnext(traj, traj, SCPP, k, i) -fnext(traj_prev, traj_prev, SCPP, k, i)))
        JtestDen = max(JtestDen, norm(fnext(traj, traj, SCPP, k, i)))
			else		# Entering new class of constraints
				JNum += JtestNum	# Add max from previous class of constraints to cost
        JDen += JtestDen
				f = fnext	# Update function handle representing current class of constraints
				JtestNum = norm(fnext(traj, traj, SCPP, k, i) - fnext(traj_prev, traj_prev, SCPP, k, i)) 	# Evaluate value of first member of new class of constraints
        JtestDen = norm(fnext(traj, traj, SCPP, k, i))
			end
		end
	end
	return JNum/JDen
end
